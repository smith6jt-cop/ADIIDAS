import ij.gui.Roi
import java.lang.Math
import qupath.lib.classifiers.pixel.PixelClassifier
import qupath.lib.images.servers.openslide.OpenslideImageServer
import qupath.lib.images.ImageData
import qupath.lib.images.servers.ImageServer
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjectTools
import qupath.lib.objects.classes.PathClass
import qupath.lib.roi.RoiTools
import qupath.lib.roi.ROIs
import qupath.lib.objects.PathObjects
import qupath.opencv.ml.pixel.PixelClassifiers
import qupath.lib.color.ColorDeconvolutionStains
import qupath.lib.images.servers.ColorTransforms
import qupath.opencv.ops.ImageOp
import qupath.opencv.ops.ImageOps
import qupath.lib.color.ColorDeconvolutionStains
import qupath.lib.color.StainVector
import qupath.lib.gui.viewer.QuPathViewer
import qupath.lib.images.servers.TransformedServerBuilder
import qupath.lib.roi.interfaces.ROI
import qupath.imagej.tools.IJTools
import qupath.lib.regions.RegionRequest
import qupath.lib.display.ChannelDisplayInfo
import qupath.opencv.ml.pixel.PixelClassifierTools
import ij.process.ImageProcessor
import qupath.lib.gui.charts.Charts
import ij.gui.HistogramWindow
import qupath.imagej.tools.PathImagePlus
import qupath.lib.gui.charts.HistogramDisplay
import qupath.lib.gui.measure.PathTableData
import ij.process.ImageStatistics
import qupath.lib.analysis.stats.Histogram
import qupath.lib.display.ImageDisplay
import java.util.concurrent.ForkJoinPool
import java.util.stream.Collectors
import java.nio.file.FileSystems
import java.nio.file.Path
import java.nio.file.Paths
import java.util.Collections

// Function to get the intensity value at a specific percentile from the histogram
def get_percentile_intensity(hist_Count, intensity_values_list, target_percentile){
    int total = 0
    for (int count : hist_Count) {
        total += count
    }

    int target = (int) (total * target_percentile)
    int cumulativeSum = 0

    for (int i = 0; i < hist_Count.size(); i++) {
        cumulativeSum += hist_Count.get(i)
        if (cumulativeSum >= target) {
            return intensity_values_list.get(i)
        }
    }
    return 0.0
}

// Function to get the percentile for a given intensity threshold
def getPercentileFromIntensity(histCount, intensityValuesList, thresholdIntensity) {
    def total = histCount.sum()
    def cumulativeSum = 0

    for (int i = 0; i < intensityValuesList.size(); i++){
        cumulativeSum += histCount[i]
        def percentile = cumulativeSum / total
        intensity = intensityValuesList.get(i)
        if (intensity >= thresholdIntensity) {
            return percentile
        }
    }
    return 1.0
}

// Function to find the first stable intensity value in the histogram, using a sliding window approach
def getFirstStableIntensity(histCount, intensityValuesList) {
    double sizeOfWindow = 10
    List previousValues = [null]*sizeOfWindow
    List averageValues = [null]*sizeOfWindow

    double averageChangeThreshold = 0.02

    int rangeIndex = 0
    int globalMax = 0
    double zeros = 0
    double previousAverage = 0
    int intensityValuesListIndex = 0
    boolean first =true
    minCount = 0
    maxCount = 0
    double sigma_max = 0
    zero_average = 0
    intensityValuesListIndex_sigma = 1
    average_average_change = 100000

    List<Integer> positiveCounts = new ArrayList<>()

    for (int i = 0; i < histCount.size(); i++) {

        int count = histCount.get(i)
        if (count > globalMax) {
            globalMax = count
        }

        if (count > 0){
            positiveCounts.add(intensityValuesList.get(i))
        }

        double intensity = intensityValuesList.get(i)

        previousValues[rangeIndex] = count
        rangeIndex = (rangeIndex + 1) % sizeOfWindow

        if (i >= sizeOfWindow-1) {
            int maxCount = previousValues.max()
            int minCount = previousValues.min()

            double currentAverage = previousValues.sum() / sizeOfWindow
            double averageChange = Math.abs(currentAverage - previousAverage) / previousAverage

            averageValues[rangeIndex] = averageChange

            if (i >= sizeOfWindow*2-1) {
                average_average_change = averageValues.sum() / sizeOfWindow
            }

            if (intensity > 0.05 && count == 0){
                zeros += 1

                if (zeros == 10 && first){
                    first = false
                    sigma_max = maxCount
                    intensityValuesListIndex = i
                }

            }

            if (average_average_change <= averageChangeThreshold &&
                    intensityValuesList.get(i) > 0.05 &&
                    !previousValues.contains(globalMax) &&
                    first) {
                first = false
                sigma_max = maxCount
                intensityValuesListIndex = i
            }
            previousAverage = currentAverage
        }
    }
    def n_biggest_index = positiveCounts.size()-5

    def nBiggestIndex = positiveCounts.size() - 5
    Double nBiggestIntensity = nBiggestIndex >= 0 ? positiveCounts[nBiggestIndex] * 10 : null

    if (nBiggestIntensity != null) {
        nBiggestIntensity = Math.round(nBiggestIntensity * 2) / 2.0;
    }
    Double stableIntensity = intensityValuesListIndex >= 0 ? intensityValuesList[intensityValuesListIndex] : null

    sigma_max = Math.log(sigma_max) / Math.log(2)
    sigma_max = Math.round(sigma_max * 2) / 2.0
    sigma_max += 7

    stableIntensity *= 10
    stableIntensity = Math.round(stableIntensity * 2) / 2.0
    stableIntensity /= 10

    return [stableIntensity, sigma_max]
}

// Function to retrieve a histogram by channel name from an image display
def Histogram getHistogramByName(String channelName, ImageDisplay imageDisplay) {
    for (ChannelDisplayInfo channel : imageDisplay.channelOptions) {
        String firstWord = channel.getName().split(" ")[0];
        if (firstWord.equals(channelName)) {
            return imageDisplay.getHistogram(channel)
        }
    }
    return null
}

// Function to calculate the ideal threshold and sigma for a specific stain channel
def ideal_thresh(server, channel) {
    // Define the stain vectors
    StainVector stain1 = StainVector.createStainVector("Glucagon", 0.79684, 0.55789, 0.23195)
    StainVector stain2 = StainVector.createStainVector("Insulin", 0.08702, 0.79214, 0.60411)
    StainVector stain3 = StainVector.createStainVector("CD3", 0.36293, 0.55989, 0.74485)

    // Create color deconvolution stains and apply them
    ColorDeconvolutionStains stains = new ColorDeconvolutionStains("3_Stains", stain1, stain2, stain3, 255, 255, 255)
    def transformedServer = new TransformedServerBuilder(server).deconvolveStains(stains).build()

    // Get the specific stain channel
    StainVector stain = stains.getStain(channel)
    String channelName = stain.getName()

    ImageData imageData = new ImageData(transformedServer)

    ImageDisplay id = new ImageDisplay()
    id.setImageData(imageData,false)
    Histogram histogram = getHistogramByName(channelName, id)

    // Calculate intensity values and the corresponding histogram bins
    double minV = histogram.getEdgeMin()
    double maxV = histogram.getEdgeMax()
    double histogramRange = histogram.getEdgeRange()
    double binWidth = histogramRange / histogram.nBins()

    List<Integer> countOfHistogram = new ArrayList<>()
    for (int histCount = 0;  histCount < histogram.nBins();histCount++){
        countOfHistogram.add(histogram.getCountsForBin(histCount))
    }

    List<Double> intensityValues = new ArrayList<>()
    for (int intensity_index = 0; intensity_index < histogram.nBins(); intensity_index++) {
        intensityValues.add(minV + binWidth * intensity_index)
    }

    // Apply log transformation to the histogram
    List<Double> logHistogram = new ArrayList<>()
    for (int i = 0; i < countOfHistogram.size(); i++) {
        if (countOfHistogram[i] > 0){
            logHistogram.add(Math.log(countOfHistogram[i]) / Math.log(2))
        } else {
            logHistogram.add(0)
        }
    }

    target_percentile = 0.996

    // Get the first stable intensity and other related parameters
    array = getFirstStableIntensity(countOfHistogram,intensityValues)
    static_thresh = array.get(0)
    target_percentile_intensity = get_percentile_intensity(countOfHistogram,intensityValues,target_percentile)
    stable_intensity = array.get(0)
    double sigma = array.get(1)
    sigma = Math.round(sigma * 2) / 2.0
    stable_intensity = Math.round(stable_intensity *20) / 20

    percentile_from_intensity = getPercentileFromIntensity(countOfHistogram,intensityValues,stable_intensity)
    print(channelName+" "+stable_intensity+" $sigma Sigma Threshold Percentile: $percentile_from_intensity")
    print("$channelName $target_percentile Percentile Intensity: $target_percentile_intensity")
    if (percentile_from_intensity < 0.99){
        print(channelName+" "+static_thresh+" Threshold Percentile Changed from: $percentile_from_intensity to $target_percentile")
        array[0] = target_percentile_intensity
        stable_intensity = array.get(0)
        stable_intensity = Math.round(stable_intensity*20) / 20

    }

    return [stable_intensity,sigma]
}

// Function to clear all annotations from the current image
def clear_annotations() {
    clearAnnotations()
}

// Function to split objects of a specific class into individual annotations
def split_class_objects(String object_class){
    resetSelection()
    selectObjects{return it.getPathClass() == getPathClass(object_class)}
    runPlugin('qupath.lib.plugins.objects.SplitAnnotationsPlugin', '{}')
}

// Function to apply a stain-based algorithm to classify and annotate image regions
def stain_algo(string_name,thresh,server, image_data,sigma) {
    if (sigma<1 ){
        sigma = 4
    } else if (sigma >16){
        sigma = 14
    }

    int stain_number = 0
    if (string_name.equals("Glucagon")){
        stain_number = 1
    }

    if (string_name.equals("Insulin")){
        stain_number = 2
    }

    print(string_name+" : "+thresh+" : "+sigma)

    // Define the stain vectors
    StainVector stain1 = StainVector.createStainVector("Glucagon", 0.79684,  0.55789,  0.23195)
    StainVector stain2 = StainVector.createStainVector("Insulin", 0.08702,  0.79214,  0.60411)
    StainVector stain3 = StainVector.createStainVector("CD3", 0.36293,  0.55989,  0.74485)

    // Create color deconvolution stains and apply them
    ColorDeconvolutionStains stains = new ColorDeconvolutionStains("3_Stains", stain1, stain2, stain3, 255, 255, 255)
    def transformedServer = new TransformedServerBuilder(server).deconvolveStains(stains, stain_number).build()
    def cal = transformedServer.getPixelCalibration()
    def resolution = cal.createScaledInstance(1, 1)
    def stain_channel = ColorTransforms.createColorDeconvolvedChannel(stains, stain_number)

    def List<ImageOp> stain_ops = new ArrayList<>()
    stain_ops.add(ImageOps.Filters.gaussianBlur(sigma))
    stain_ops.add(ImageOps.Threshold.threshold(thresh))

    Map<Integer, PathClass> stain_classifications = new LinkedHashMap<>()
    stain_classifications.put(0, null)
    stain_classifications.put(1, getPathClass(string_name))

    // Create a classifier for the stain and apply it
    def stain_op = ImageOps.Core.sequential(stain_ops)
    def stain_transformer = ImageOps.buildImageDataOp(stain_channel).appendOps(stain_op)
    stain_classifier = PixelClassifiers.createClassifier(
            stain_transformer,
            resolution,
            stain_classifications
    )

    print("Running Created Classifier")
    PixelClassifierTools.createAnnotationsFromPixelClassifier(image_data, stain_classifier, 20.0, 0.0)

    print "Finished Splitting $string_name objects"
    resetSelection()
}

// Function to check if an annotation's centroid is within certain bounds
def checkoverlap_centroid(centroidROI, bounds) {
    x = centroidROI.getCentroidX()
    y = centroidROI.getCentroidY()
    def (x_min, y_min, x_max, y_max) = bounds

    return x_min <= x && x < x_max && y_min <= y && y < y_max
}

// Function to merge annotations that overlap within specific bounds
def mergeAnnotationsWithinBounds() {
    int index_1 = 0
    loopList = getObjects() {it.getPathClass() == getPathClass("Islet")}
    while (index_1 < loopList.size()) {
        loopList = getObjects() {it.getPathClass() == getPathClass("Islet")}
        loopsize = loopList.size()

        if (index_1 >= loopList.size()) {
            break
        }

        def initial = loopList[index_1]
        def mergedAnnotations = [initial]
        def initialROI = initial.getROI()
        x_min = initialROI.getBoundsX()
        y_min = initialROI.getBoundsY()
        x_max = initialROI.getBoundsX() + initialROI.getBoundsWidth()
        y_max = initialROI.getBoundsY() + initialROI.getBoundsHeight()

        targetAnnotationList = getObjects() {
            it.getPathClass() == getPathClass("Islet") && (checkoverlap_centroid(it.getROI(), [x_min, y_min, x_max, y_max]))
                    && !it.equals(initial)
        }

        for (target : targetAnnotationList) {
            mergedAnnotations.add(target)
        }

        if (mergedAnnotations.size() > 1) {
            def ins_area = 0
            def glu_area = 0
            for (measurement_annotation : mergedAnnotations) {
                ins_area += measurement_annotation.getMeasurementList().get("Insulin")
                glu_area += measurement_annotation.getMeasurementList().get("Glucagon")
            }
            selectObjects(mergedAnnotations)
            mergeSelectedAnnotations()
            def new_merged = getSelectedObject()
            def new_merged_list = new_merged.getMeasurementList()
            new_merged_list.put("Insulin", ins_area)
            new_merged_list.put("Glucagon", glu_area)
            if (loopsize > loopList.size()) {
                index_1 -= (loopsize-loopList.size())
                if (index_1 < 0) {
                    index_1 = 0
                }
            }
        }
        else {
            index_1++
        }
    }
}

// Function to delete all annotations of a specific class
def delete_annotation(name){
    selectObjects{return it.getPathClass() == getPathClass(name)}
    clearSelectedObjects(false)
}

// Function to check if a class has more than a certain number of annotations, and delete them if so
def listBiggerThanLimit(name,count_limit){
    list = getObjects() {return it.getPathClass() == getPathClass(name)}
    if (list.size() > count_limit){
        selectObjects(list)
        clearSelectedObjects(false)
        return true
    } else {
        return false
    }
}

// Function to delete small annotations below a certain size threshold
def delete_small(name, size, scale_param){
    selectObjects{return it.getPathClass() == getPathClass(name) && it.getROI().getScaledArea(scale_param,scale_param) < size}
    clearSelectedObjects(false)
}

// Function to calculate and store hole measurements (e.g., hole count and hole area ratio) for each annotation
def makeHoleMeasurements(name){
    loopList = getObjects() { it.getPathClass() == getPathClass(name) }
    for (annotation : loopList) {
        measurements = annotation.getMeasurementList()
        annotation_area = annotation.getROI().getArea()
        polygons = RoiTools.splitAreaToPolygons(annotation.getROI())

        hole_count = polygons[0].size()
        holes_area = 0
        for (hole: polygons[0]) {
            holes_area += hole.getArea()
        }

        hole_ratio = holes_area/annotation_area

        measurements.put("Holes", hole_count)
        measurements.put("Hole Ratio", hole_ratio)
    }
}

// Function to calculate and store circularity measurements for each annotation, returning the mean circularity for annotations larger than a cutoff
def makeCircularityMeasurements(name, cutoff){
    loopList = getObjects() { it.getPathClass() == getPathClass(name) }
    double mean = 0
    double count = 0
    for (annotation : loopList) {
        measurements = annotation.getMeasurementList()
        annotation_area = annotation.getROI().getArea()
        double perimeter = annotation.getROI().getLength()

        if (!annotation_area.isNaN() && !perimeter.isNaN() && perimeter != 0) {
            circularity = 4 * Math.PI * (annotation_area / (perimeter ** 2))
            measurements.put("Circularity", circularity)
        } else {
            circularity = 0
            measurements.put("Circularity", circularity)
        }

        if (annotation_area >= cutoff){
            mean += circularity
            count += 1
        }
    }
    if (count>0 && mean > 0){
        mean = mean/count
        return mean
    } else {
        return 0
    }
}

// Function to calculate and store diameter measurements for each annotation, returning the mean diameter for annotations larger than a cutoff
def makeDiameterMeasurements(name, cutoff){
    loopList = getObjects() { it.getPathClass() == getPathClass(name) }
    double mean = 0
    double count = 0
    for (annotation : loopList) {
        measurements = annotation.getMeasurementList()
        annotation_area = annotation.getROI().getArea()
        double diameter = 2*Math.sqrt(annotation_area/Math.PI)
        if (annotation_area >= cutoff){
            mean += diameter
            count +=1
        }

        measurements.put("Diameter", diameter)
    }
    if (count>0 && mean > 0){
        mean = mean/count
        return mean
    } else {
        return 0
    }
}

// Function to check if two bounding boxes overlap
def checkOverlap(box1, box2) {
    def (x1_min, y1_min, x1_max, y1_max) = box1
    def (x2_min, y2_min, x2_max, y2_max) = box2

    if (x1_max < x2_min || x2_max < x1_min) {
        return false  // No horizontal overlap
    }

    if (y1_max < y2_min || y2_max < y1_min) {
        return false  // No vertical overlap
    }

    return true  // Overlap detected
}

// Function to calculate the total area of a specific marker class
def marker_area_collection_function(string_name) {
    double area = 0.0
    string_List = getObjects() {return it.getPathClass() == getPathClass(string_name)}
    for(object:string_List){
        annotation_area = object.getROI().getArea()
        if (!(annotation_area.isNaN())) {
            area += annotation_area
        }
    }

    return area > 0 ? area : 0.0
}

// Function to calculate the total area of insulin and glucagon for islets
def area_collection_function(string_name) {
    double area = 0.0
    IsletList = getObjects() { it.getPathClass() == getPathClass("Islet")}
    for(object:IsletList){
        area += object.getMeasurementList().get(string_name)
    }

    return area > 0 ? area : 0.0
}

// Function to assign the percentage of the bounding box area that the annotation occupies
def assign_bounds_percent(string_name) {
    loopList = getObjects() { it.getPathClass() == getPathClass(string_name) }
    if (!(loopList.isEmpty())) {
        for (annotation : loopList) {
            measurements = annotation.getMeasurementList()
            annotation_area = annotation.getROI().getArea()
            width = annotation.getROI().getBoundsWidth()
            height = annotation.getROI().getBoundsHeight()
            bounds_area = width * height

            if (!bounds_area.isNaN() && !annotation_area.isNaN()) {
                area_percent = 100 * (annotation_area / bounds_area)
                measurements.put("Bounds Area Percent", area_percent)
            } else {
                measurements.put("Bounds Area Percent", 0)
            }
        }
    }
}

// Function to calculate and store the insulin to glucagon ratio for each islet
def makeInsulinToGlucagonRatios(name){
    loopList = getObjects() { it.getPathClass() == getPathClass(name) }

    for (annotation : loopList) {
        measurements = annotation.getMeasurementList()
        insArea = measurements.get("Insulin")
        gluArea = measurements.get("Glucagon")

        if (insArea==0 || gluArea==0){
            insGluRatio = 0
        } else {
            insGluRatio = insArea / gluArea
        }
        measurements.put("Insulin/Glucagon Ratio", insGluRatio)
    }
}

// Function to convert insulin and glucagon annotations to islet annotations
def convert_endocrine_to_islet(scale){
    def insulinList = getObjects(){it.getPathClass() == getPathClass("Insulin") }
    def glucagonList = getObjects(){ it.getPathClass() == getPathClass("Glucagon") }
    def imageData = getCurrentImageData()

    int processors = Runtime.getRuntime().availableProcessors()
    ForkJoinPool customThreadPool = new ForkJoinPool(processors)
    def newInsulinObjects = []
    def newGlucagonObjects = []

    try {
        customThreadPool.submit({
            newInsulinObjects = insulinList.parallelStream().map { insulin ->
                def pathObject = PathObjects.createAnnotationObject(insulin.getROI(), getPathClass("Islet"))
                def measurements = pathObject.getMeasurementList()
                measurements.put("Insulin", insulin.getROI().getScaledArea(scale, scale))
                measurements.put("Glucagon", 0)
                return pathObject
            }.collect(Collectors.toList())

            newGlucagonObjects = glucagonList.parallelStream().map { glucagon ->
                def pathObject = PathObjects.createAnnotationObject(glucagon.getROI(), getPathClass("Islet"))
                def measurements = pathObject.getMeasurementList()
                measurements.put("Insulin", 0)
                measurements.put("Glucagon", glucagon.getROI().getScaledArea(scale, scale))
                return pathObject
            }.collect(Collectors.toList())
        }).get()
    } finally {
        addObjects(newInsulinObjects)
        addObjects(newGlucagonObjects)
        removeObjects(insulinList,false)
        removeObjects(glucagonList,false)
        customThreadPool.shutdown()
    }
}

// Function to group annotations of islets that are large enough
def group_algo (string_name,scale_param, size_param) {
    loopList = getObjects() {
        it.isAnnotation() &&
                it.getROI().getScaledArea(scale_param,scale_param) >= size_param
    }

    if (loopList.isEmpty()) {
        print "No annotation object found for grouping."
    } else if (loopList.size() < 2) {
        print "Not enough $string_name annotation objects found for grouping."
    } else {
        mergedAnnotations = []
        int boxDistance = 300 / scale_param
        int boundsDistance = 0
        int index_1 = 0

        while (index_1 < loopList.size()) {

            loopList = getObjects() {
                it.isAnnotation() &&
                        it.getROI().getScaledArea(scale_param,scale_param) >= size_param}
            loopsize = loopList.size()

            if (index_1 >= loopList.size()) {
                break
            }

            def initial = loopList[index_1]
            def mergedAnnotations = [initial]
            def initialROI = initial.getROI()
            x = initialROI.getCentroidX()
            y = initialROI.getCentroidY()

            comparisonList = getObjects() {
                it.isAnnotation() &&
                        Math.abs(it.getROI().getCentroidX() - x) <= boxDistance &&
                        Math.abs(it.getROI().getCentroidY() - y) <= boxDistance &&
                        !it.equals(initial)
            }

            for (comparison : comparisonList) {
                if (RoiTools.getBoundaryDistance(initialROI, comparison.getROI(), scale_param, scale_param) <= boundsDistance) {
                    mergedAnnotations.add(comparison)
                }
            }

            if (mergedAnnotations.size() == 1) {
                index_1++
            } else {
                ins_area = 0
                glu_area = 0
                for (measurement_annotation : mergedAnnotations){
                    ins_area += measurement_annotation.getMeasurementList().get("Insulin")
                    glu_area += measurement_annotation.getMeasurementList().get("Glucagon")
                }
                selectObjects(mergedAnnotations)
                mergeSelectedAnnotations()
                new_merged = getSelectedObject()
                new_merged_list = new_merged.getMeasurementList()
                new_merged_list.put("Insulin", ins_area)
                new_merged_list.put("Glucagon", glu_area)
                if (loopsize > loopList.size()) {
                    index_1 -= (loopsize-loopList.size())
                    if (index_1 < 0) {
                        index_1 = 0
                    }
                }
            }
        }
    }
}

// Function to draw lines between annotations of a specific class and calculate their distances
def line_algo(loopList,server,scale, string_name, size_cutoff){
    linesDrawn = []
    int clusters = 0

    lineLengths = []

    if (string_name.equals("Islets")){
        lineType = "Islets Line"
    }

    if (loopList.isEmpty() || loopList.size() < 2) {
        if (string_name=="Glucagon"){
            print "Insufficient Insulin Negative annotation objects found for line algo."
        } else {
            print "Insufficient $string_name annotation objects found for line algo."
        }
    } else {
        for (def annotation : loopList) {
            annotationRoi = annotation.getROI()
            nearestAnnotation = null
            double minDistance = Double.MAX_VALUE

            x = annotationRoi.getCentroidX()
            y = annotationRoi.getCentroidY()

            int line_box_distance = 3000 / scale

            comparisonList = getObjects() {
                it.isAnnotation() &&
                        Math.abs(it.getROI().getCentroidX() - x) <= line_box_distance &&
                        Math.abs(it.getROI().getCentroidY() - y) <= line_box_distance &&
                        loopList.contains(it) &&
                        it.getROI().getScaledArea(scale,scale) >= size_cutoff &&
                        !it.equals(annotation)
            }

            for (def otherAnnotation : comparisonList) {
                otherAnnotationRoi = otherAnnotation.getROI()
                distance = RoiTools.getCentroidDistance(annotationRoi, otherAnnotationRoi, 1.0, 1.0)
                line = ROIs.createLineROI(annotationRoi.getCentroidX(), annotationRoi.getCentroidY(), otherAnnotationRoi.getCentroidX(), otherAnnotationRoi.getCentroidY(), annotationRoi.getImagePlane())
                if (distance < minDistance) {
                    minDistance = distance
                    nearestAnnotation = otherAnnotation
                }
            }
            if (nearestAnnotation != null){
                if (!(linesDrawn.contains([annotationRoi, nearestAnnotation.getROI()]) || linesDrawn.contains([nearestAnnotation.getROI(), annotationRoi]))) {
                    line = ROIs.createLineROI(annotationRoi.getCentroidX(), annotationRoi.getCentroidY(), nearestAnnotation.getROI().getCentroidX(), nearestAnnotation.getROI().getCentroidY(), annotationRoi.getImagePlane())

                    pathClass = getPathClass(lineType)

                    lineLength = line.getScaledLength(scale, scale)
                    if (lineLength <10000) {
                        lineLengths.add(lineLength)
                        pathObject = PathObjects.createAnnotationObject(line,pathClass)
                        imageData = getCurrentImageData()
                        imageData.getHierarchy().addObject(pathObject)

                        linesDrawn << [annotationRoi, nearestAnnotation.getROI()]
                    }
                } else {
                    clusters += 1
                }
            }
        }

        String name = server.getMetadata().getName()
        name = name.toString()
        name = name.split("\\.", 0)
        name = name[0]

        int lineCount = lineLengths.size()
        double meanLength = 0.0
        double median = 0

        if (lineCount > 0) {
            meanLength = lineLengths.stream().mapToDouble(Double::doubleValue).average().orElse(0)
            List<Double> sortedLengths = new ArrayList<>(lineLengths)
            Collections.sort(sortedLengths)

            int index = lineCount / 2

            if (lineCount == 1) {
                median = sortedLengths.get(0)
            } else {
                if (lineCount % 2 != 0 && lineCount != 1) {
                    median = (sortedLengths.get(index) + sortedLengths.get(index + 1)) / 2
                } else {
                    median = sortedLengths.get(index)
                }
            }
        }
        return new double[]{clusters, lineCount, meanLength, median}
    }
    return new double[]{0, 0, 0, 0}
}

// Main function that orchestrates the analysis process
def main() {
    // Initialize necessary parameters for environment
    setImageType('BRIGHTFIELD_OTHER')
    active_server = getCurrentServer()
    def scale = 0
    if (active_server instanceof OpenslideImageServer) {
        OpenslideImageServer openSlideServer = (OpenslideImageServer) active_server
        Map<String, String> properties = openSlideServer.getProperties()
        def metadata = properties['originalMetadata']
        scale = metadata['pixelWidthMicrons']
        print("Pixel Size Set: "+scale)
        setPixelSizeMicrons(scale, scale)
    }

    ImageData activeImageData = getCurrentImageData()

    String name = active_server.getMetadata().getName().toString().split("\\.", 0)[0]
    print(name)

    // Delete any previous annotations
    clear_annotations()

    // Gets the predefined tissue classifier from the Qupath
    print("Tissue Classifier")
    createAnnotationsFromPixelClassifier("Tissue_Classifier_1", 0, 0.0)
    tissueAnnotations = getObjects() { it.getPathClass() == getPathClass("Tissue") }

    //Collecting the Area from the annotations qupath recognised.
    print("Collecting Tissue Area")
    double tissue_area = tissueAnnotations.isEmpty() ? 0 : tissueAnnotations.collect { it.getROI().getScaledArea(scale, scale) }.sum()
    //delete_annotation("Tissue")

    glu_thresh = ideal_thresh(active_server, 1)
    ins_thresh = ideal_thresh(active_server, 2)
    print("Glucagon Ideal Thresh: " + glu_thresh.get(0))
    print("Insulin Ideal Thresh: " + ins_thresh.get(0))
    print("Glucagon Ideal Sigma: " + glu_thresh.get(1))
    print("Insulin Ideal Sigma: " + ins_thresh.get(1))

    def insulinLimit = true
    def insulinEndocrineAreaPercent = 0
    ins_thresh_static = ins_thresh.get(0)
    if (ins_thresh_static < 0.05) {
        ins_thresh_static = 0.3
    } else if (ins_thresh_static >= 0.4) {
        ins_thresh_static = 0.2
    }
    while (insulinLimit || insulinEndocrineAreaPercent > 10) {
        delete_annotation("Insulin")
        print("Insulin Deconvolution")
        stain_algo("Insulin", ins_thresh_static, active_server, activeImageData, ins_thresh.get(1))
        insulinLimit = listBiggerThanLimit("Insulin", 12000)
        insulinEndocrineAreaPercent = marker_area_collection_function("Insulin") * 100 / tissue_area
        print("Insulin Endocrine Area Percent: " + insulinEndocrineAreaPercent)
        if (insulinEndocrineAreaPercent > 100) {
            ins_thresh_static += 0.05
        }
        ins_thresh_static += 0.05
    }

    def glucagonLimit = true
    def glucagonEndocrineAreaPercent = 0
    glu_thresh_static = glu_thresh.get(0)
    if (glu_thresh_static < 0.05) {
        glu_thresh_static = 0.3
    } else if (glu_thresh_static >= 0.4) {
        glu_thresh_static = 0.2
    }
    while (glucagonLimit || glucagonEndocrineAreaPercent > 10) {
        delete_annotation("Glucagon")
        print("Glucagon Deconvolution")
        stain_algo("Glucagon", glu_thresh_static, active_server, activeImageData, glu_thresh.get(1))
        glucagonLimit = listBiggerThanLimit("Glucagon", 12000)
        glucagonEndocrineAreaPercent = marker_area_collection_function("Glucagon") * 100 / tissue_area
        print("Glucagon Endocrine Area Percent: " + glucagonEndocrineAreaPercent)
        if (glucagonEndocrineAreaPercent > 100) {
            glu_thresh_static += 0.05
        }
        glu_thresh_static += 0.05
    }

    split_class_objects("Glucagon")
    split_class_objects("Insulin")

    print("Making Conversions")
    convert_endocrine_to_islet(scale)

    // Group regions large enough to create full sized islets
    print("Grouping")
    group_algo("Islet", scale, 0)
    IsletList = getObjects() { it.getPathClass() == getPathClass("Islet") }

    mergeAnnotationsWithinBounds()

    double cutoff = 1000

    // Create lists for measuring the endocrine area
    print("Area Collection Functions")
    resetSelection()

    // Initiate Area Collection
    double glucagon_area = area_collection_function("Glucagon")
    double insulin_area = area_collection_function("Insulin")
    double total_endocrine_area = glucagon_area + insulin_area

    // Assign Formulaic Measurements
    assign_bounds_percent("Islet")
    makeHoleMeasurements("Islet")
    meanDiameter = makeDiameterMeasurements("Islet", cutoff)
    meanCircularity = makeCircularityMeasurements("Islet", cutoff)
    makeInsulinToGlucagonRatios("Islet")

    // Retrieve true INS+ Islets and INS- Islets
    print("Retreiving Annotations")

    normal_IsletList = getObjects() {
        it.getPathClass() == getPathClass("Islet") &&
                it.getROI().getScaledArea(scale, scale) >= cutoff
    }

    insulinPositiveIsletList = getObjects() {
        it.getPathClass() == getPathClass("Islet") &&
                it.getROI().getScaledArea(scale, scale) >= cutoff &&
                it.getMeasurementList().get("Insulin") > 0
    }

    insulinNegativeIsletList = getObjects() {
        it.getPathClass() == getPathClass("Islet") &&
                it.getROI().getScaledArea(scale, scale) >= cutoff &&
                it.getMeasurementList().get("Insulin") == 0 &&
                it.getMeasurementList().get("Glucagon") > 0
    }
    IsletCount = normal_IsletList.size()
    insulinPositiveIsletCount = insulinPositiveIsletList.size()
    insulinNegativeIsletCount = insulinNegativeIsletList.size()

    // Retrieve Islet-Hormone Cluster Marker Combinations
    insulinPositiveGlucagonPositiveIsletList = getObjects() {
        it.getPathClass() == getPathClass("Islet") &&
                it.getROI().getScaledArea(scale, scale) >= cutoff &&
                it.getMeasurementList().get("Insulin") > 0 &&
                it.getMeasurementList().get("Glucagon") > 0
    }

    insulinNegativeGlucagonPositiveIsletList = getObjects() {
        it.getPathClass() == getPathClass("Islet") &&
                it.getROI().getScaledArea(scale, scale) >= cutoff &&
                it.getMeasurementList().get("Insulin") == 0 &&
                it.getMeasurementList().get("Glucagon") > 0
    }

    insulinPositiveGlucagonNegativeIsletList = getObjects() {
        it.getPathClass() == getPathClass("Islet") &&
                it.getROI().getScaledArea(scale, scale) >= cutoff &&
                it.getMeasurementList().get("Insulin") > 0 &&
                it.getMeasurementList().get("Glucagon") == 0
    }

    insulinPositiveGlucagonPositiveHormoneClusterList = getObjects() {
        it.getPathClass() == getPathClass("Islet") &&
                it.getROI().getScaledArea(scale, scale) < cutoff &&
                it.getMeasurementList().get("Insulin") > 0 &&
                it.getMeasurementList().get("Glucagon") > 0
    }

    insulinNegativeGlucagonPositiveHormoneClusterList = getObjects() {
        it.getPathClass() == getPathClass("Islet") &&
                it.getROI().getScaledArea(scale, scale) < cutoff &&
                it.getMeasurementList().get("Insulin") == 0 &&
                it.getMeasurementList().get("Glucagon") > 0
    }

    insulinPositiveGlucagonNegativeHormoneClusterList = getObjects() {
        it.getPathClass() == getPathClass("Islet") &&
                it.getROI().getScaledArea(scale, scale) < cutoff &&
                it.getMeasurementList().get("Insulin") > 0 &&
                it.getMeasurementList().get("Glucagon") == 0
    }

    print("Line Algo Start")
    double[] allVariables = line_algo(normal_IsletList, active_server, scale,"Islets",cutoff)
    int Islet_Clusters = allVariables[0]
    int Islet_LineCount = allVariables[1]
    double Islet_MeanLength = allVariables[2]
    double Islet_Median = allVariables[3]

    endocrine_to_tissue_ratio = total_endocrine_area * 100 / tissue_area

    double islet_sum = insulinPositiveGlucagonPositiveIsletList.size() +
            insulinNegativeGlucagonPositiveIsletList.size() +
            insulinPositiveGlucagonNegativeIsletList.size()

    islet_mm_density = (islet_sum * 1000000 / tissue_area)
    def outputDir = buildFilePath(PROJECT_BASE_DIR, 'export')
    mkdirs(outputDir)

    def path = buildFilePath(outputDir, "Output_FM_Beauty.txt")

    boolean writeHeader = !new File(path).exists()
    try (FileWriter writer = new FileWriter(path, true)) {
        if (writeHeader) {
            writer.write("ImageID" + "," +
                    "tissue_area_µm^2" + "," + "endocrine_area_total" + "," +
                    "endocrine_area_%_tissue_area" + "," + "islet_density_per 1mm^2 tissue" + "," +
                    "mean_islet_diameter" + "," + "mean_islet_circularity" + "," +
                    "ins+_islets" + "," +
                    "ins-_islets" + "," +
                    "ins+/gcg+_clusters" + "," + "ins-/gcg+_clusters" + "," +
                    "ins+/gcg-_clusters" + "," + "ins+/gcg+_islets" + "," +
                    "ins-/gcg+_islets" + "," + "ins+/gcg-_islets" + "," +
                    "Interislet_distance_Mean" + "," +
                    "Interislet_distance_Median" + "," + "Insulin_area" + "," + "Glucagon_area")
        }
        writer.write("\n" + name + "," +
                tissue_area + "," + total_endocrine_area + "," + endocrine_to_tissue_ratio + "," +
                islet_mm_density + "," + meanDiameter + "," + meanCircularity + "," +
                insulinPositiveIsletCount + "," +
                insulinNegativeIsletCount + "," +
                insulinPositiveGlucagonPositiveHormoneClusterList.size() + "," +
                insulinNegativeGlucagonPositiveHormoneClusterList.size() + "," +
                insulinPositiveGlucagonNegativeHormoneClusterList.size() + "," +
                insulinPositiveGlucagonPositiveIsletList.size() + "," +
                insulinNegativeGlucagonPositiveIsletList.size() + "," +
                insulinPositiveGlucagonNegativeIsletList.size() + "," +
                Islet_MeanLength + "," + Islet_Median + "," + insulin_area + "," + glucagon_area)
    }
    print("Execution Successful")
}
main()
