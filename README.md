# ADIIDAS (Algorithm Determining Inter-Islet Distance And Size)

## This code is used to perform semi-automated image analysis of whole-slide images of pancreas tissue sections stained with classic immunohistochemistry methods against insulin in red and glucagon in blue. Images in SVS format are imported in QuPath. The code can be executed for individual images or as a batch to the entire QuPath project to obtain measurements regarding insulin and glucagon areas, annotations counts and sizes as well as distance to the closest annotation (≥1000µm2 area).

### Critical parameters for this dataset:

&emsp;- Tissue_classifier_1 is predefined, as attached, with a sigma values of 2 and threshold of 230

&emsp;- StainVector stain1 = StainVector.createStainVector("Glucagon", 0.79684, 0.55789, 0.23195)

&emsp;- StainVector stain2 = StainVector.createStainVector("Insulin", 0.08702, 0.79214, 0.60411)

&emsp;- StainVector stain3 = StainVector.createStainVector("CD3", 0.36293, 0.55989, 0.74485)


The code starts with setting up the environment, setting the image type, extracting the metadata and pixel size, and clearing any previous annotations if applicable. The predefined tissue classifier is used to detect the total tissue on the slide and then the total tissue area is measured.

### Main functions:

- Ideal thresholds for insulin and glucagon are calculated using histograms of pixel intensity
  - Histogram of pixel intensity is created
    - Ideal thresholds are selected using the sliding window approach (getsfirststableintensity function), which determines the point in the histogram where the intensity values start to stabilize
    - The threshold percentile is checked to match a target value, otherwise adjusted to match the target percentiles

&emsp;- Stained regions are classified and annotated as insulin and glucagon using thresholds identified in the previous step

        - Stain vectors are defined based on OD values for red and blue, as specified in the critical parameter section
        - Color deconvolution is applied to separate the red and blue stains from background/other stains
        - Gaussian blur is applied based on sigma and threshold values identified using the histograms of pixel intensity in order to smoothen the image (reduce noise) and separate stained regions from background
        - Pixel classifiers are created and run
        - Annotations are created wherever the thresholds defined above are being met
        - Annotations are refined and merged when stains are in contact using the “group_algo” function

&emsp;- Distance between islets is calculated using the “line_algo” function

        - Centroids for each islets are created and x and y coordinates for each annotation are determined
        - A line annotation is created between the centroid of an annotation and the nearest annotations, and the length of the line is measured to determine the distance to the closest islet

&emsp;- Measurements are exported.
