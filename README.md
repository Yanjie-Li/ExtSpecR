
# ExtSpecR
![Screen](/images/22.png)

<!-- badges: start -->
<!-- badges: end -->

## Citation

If you use this project, please cite the following paper
### BibTeX

```
@article{
doi:10.34133/plantphenomics.0103,
author = {Zhuo Liu  and Mahmoud Al-Sarayreh  and Cong Xu  and Federico Tomasetto  and Yanjie Li },
title = {<i>ExtSpecR</i>: An R Package and Tool for Extracting Tree Spectra from UAV-Based Remote Sensing},
journal = {Plant Phenomics},
volume = {5},
number = {},
pages = {0103},
year = {2023},
doi = {10.34133/plantphenomics.0103},
URL = {https://spj.science.org/doi/abs/10.34133/plantphenomics.0103},
eprint = {https://spj.science.org/doi/pdf/10.34133/plantphenomics.0103}}
```
### Infographic

![Infographic](/images/D-23-00146.jpg)

## Video Tutorial
#### Check out our video tutorial for a quick introduction:
 
 
[![Watch the video](/images/1679059768729.png)](https://youtu.be/NVB1kq8AsK4)

 

## Preface
Welcome to our ExtSpecR R package! This is a shiny application specifically designed to extract spectral information from a single tree, and I also provide a sample to show what the data looks like after extraction. First of all, you need point cloud data with precise positioning information, which is used to segment each individual plant in a large forest area. You also need multispectral or hyperspectral information of that forest area.


Before installing the package, you need to make sure that you have installed  `rtools`,You can check and install `rtools` using the provided link:
## [Download rtools](https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html).  

We provided two ways to use ExtSpecR pacakge, first is the normal way:
## (First way) Install package
Before installing my package, please make sure you have installed all the required packages that my package depends on. The following R packages are needed:
`Biobase`, `RStoolbox`, `EBImage`,`shinythemes`,`shinyWidgets`,`terra`,`shinyjs`,`RCSF`,`DT`,`shinydashboard`,`stars`,`colorspace`,`readr`,`sfheaders`,`sf`,`exactextractr`,`lidR`,`tidyverse`,`viridis`,`rgdal`,`tictoc`,`ggrepel`,`raster`,`tools`,`rasterVis`,`data.table`,
If you do not have any of these packages installed, you can use the following code to check and install them:
``` r
 # First, check and install required packages
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("bleutner/RStoolbox")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

pkgs <- c("Biobase", "EBImage")
if (length(setdiff(pkgs, rownames(installed.packages()))) > 0)
  BiocManager::install(pkgs, ask = FALSE, update = FALSE)

packages <- c("shinythemes", 'shinyWidgets', 'terra', 'shinyjs', 'RCSF', 'DT',
              "shinydashboard", 'stars', 'colorspace', 'readr',
              'sfheaders', 'sf', 'exactextractr', 'lidR', 'tidyverse',
              'viridis', 'rgdal', 'tictoc', 'ggrepel',
              'raster', 'tools', 'rasterVis', 'data.table')

new.packages <- packages[!(packages %in% utils::installed.packages()[,"Package"])]
if(length(new.packages)) utils::install.packages(new.packages, quiet = F)

# Once all dependencies are installed, install ExtSpecR
devtools::install_github("Yanjie-Li/ExtSpecR")


```
 

## Example

This is a basic example:

``` r
library(ExtSpecR)
## basic example code
ExtSpecR_app()

```
## (Second way) The R-portable
If you encounter any issues during installation, you can try the following method: 
 
### 1 Download the R-portable compressed file from the link provided: 
### [Download R-portable package](https://ln5.sync.com/dl/1d8587200/aubeg7ib-x7ia9bx5-r86f7qqy-2rvuygaf)
### 2 Unzip the package to a desired location.
### 3 Navigate to the extracted directory and double-click on "R-Portable.exe" to open the R environment.
### 4 In the R console, load the ExtSpecR package using the following command:

``` r
library(ExtSpecR)
```
### 5 You can then use the functions in ExtSpecR package, for example:
 
``` r
ExtSpecR_app()
```
This will start the ExtSpecR Shiny application.
Note that the advantage of using the R-portable package is that it comes with all the necessary packages preinstalled, so you don't have to install any additional packages or worry about conflicting package versions. Also, using R-portable does not interfere with your existing system environment.

I hope you can use my R package without any problems! If you encounter any problems, please do not hesitate to contact me and I will do my best to help you.
 
## brief introduction

Let me give you a quick overview of our ExtSpecR_app() website ![Screen](/images/fiugre1.PNG) We have three main sections: Introduction, VIs and Examples, and Tree Phenotyping.

In the Introduction section, we've included a handy flowchart to guide you through the steps of using our application. It's a great way to get a feel for what to expect and how to get started.

Moving on to VIs and Examples, we have two sub-menus: Generate VIs and Example Data. The first allows you to upload your own five-band TIF raster images from multispectral sensors and calculate different vegetation indices while visualising them on the image display page. Just select the index you want to see and go to the image viewer page, and you'll get a false-colour image and corresponding vegetation index map. Be patient, however, as the data may take some time to process due to the size of the TIF files. When it's done, you can download the resulting images with the click of a button.


The workflow for the VIs calculation is shown in the flowing graph:
![Screenshot](/images/VIs.png)

The page is showing like this:
![Screenshot](/images/figurevis.PNG)

In the example data section, we've included some data from my own research, extracted with my ExtSpecR package, which includes tree heights, crown area, and crown layer vegetation indices from research site. You can explore each tree's data by viewing the images and checking out the data info page.

 
Next is the most important part of this app, tree phenotyping. It includes two sub-menus: Data process and Download. The main focus is on the Data process menu, which includes four pages: Upload and Draw ROI, Segmentation, LAS plot information, and Extraction and Visualization. ![Screenshot](/images/treephno.png)

On the Upload and Draw ROI page, users need to upload point cloud data and corresponding raster images data. After uploading, click on the Draw ROI polygon button, and a separate R graphic page will appear. Users can select the area of interest by clicking with the left mouse button and finishing by clicking the right mouse button. The uploaded point cloud data, ROI-selected point cloud data, and image spectral information will be displayed on the Point cloud information and Spectral information sub-page.

After selecting the ROI, users can optimize the Segmentation process by adjusting the hmin and ws information on the Segmentation page![Screenshot](/images/seg.png). After selecting the settings, click on the Start segmentation button, and a small window will appear in the lower right corner of the page, indicating that the calculation is in progress. When the window disappears, two images will appear, a false-color image and the polygon of each tree extracted from the selected ROI area. Users can adjust the TreeID number to locate different tree positions.

After the Segmentation process is complete, users can use the LAS plot information page![Screenshot](/images/las.png) to plot the point cloud information, including the original point cloud data, ROI point cloud data, and segmented trees. The plotted single tree is the one with the TreeID number selected during Segmentation. Users can adjust the point size slider to plot information about different-sized point clouds.

The final step is Extraction and Visualization![Screenshot](/images/exdt.png). After the Segmentation process is complete, users can click on the Extract data button on the Extraction and Visualization page to extract image spectral data on a tree-by-tree basis. Once extracted, the page displays the different data visualizations of each tree's crown, including Z tree height, area crown width, and image spectral information. Users can adjust the Height limits slider to filter the crown spectral information under different height conditions.

Finally, after completing all the steps, users can click on the Download sub-menu to download the final data![Screenshot](/images/down.png), which is available as ".rds" data that can be read in R language. Users can also download the polygon data of each tree extracted during Segmentation. The right-side page shows the data format in an intuitive way.


## data

please download the example data from the bellowed link:[example data](https://ln5.sync.com/dl/d6899c6f0/3g32725x-b85yuvm3-ba68kfre-jewun6fk)

