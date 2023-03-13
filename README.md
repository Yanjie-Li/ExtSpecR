
# ExtSpecR
![Screen](/images/22.png)

<!-- badges: start -->
<!-- badges: end -->

Welcome to my ExtSpecR R package! The  is a shiny app for single tree multi-spectral extraction This is a shiny app specially used to extract spectral information of a single tree, and I also provide a sample to show what the data looks like after extraction.
First of all, you need point cloud data with precise positioning information, which is used to segment each individual plant in a large area of forest land. In addition, you need to have multi-spectral or hyperspectral information of this forest land.

There are two ways to use ExtSpecR packages, first is the normal way:
## (First way) Install package
Before you install my package, please make sure you have installed some required packages that my package depends on. The following R packages are required:
`Biobase`,`EBImage`,`shinythemes`,`shinyWidgets`,`terra`,`shinyjs`,`RCSF`,`DT`,`shinydashboard`,`stars`,`colorspace`,`readr`,`sfheaders`,`sf`,`exactextractr`,`lidR`,`tidyverse`,`viridis`,`rgdal`,`tictoc`,`ggrepel`,`raster`,`tools`,`rasterVis`,`data.table`,
If you have not installed any of these packages, you can use the following code to check and install them:
``` r
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
if(length(new.packages)) utils::install.packages(new.packages, quiet = TRUE)

```
Once you have installed all the required packages, you can proceed with the installation of the development version of ExtSpecR  package from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Yanjie-Li/ExtSpecR")
```
Enjoy using my R package!
Note:`rtools` required, `rtools` should install manually from: [rtools](https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html).  

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(ExtSpecR)
## basic example code
ExtSpecR_app()

```
## (Second way) The R-portable
If you encounter any issues during installation, you can try the following method: download the R-portable compressed file[R-portable package](https://ln5.sync.com/dl/1d8587200/aubeg7ib-x7ia9bx5-r86f7qqy-2rvuygaf), which includes all the required packages, including my ExtSpecR package. After downloading, you can unzip the package, and then click on the "R-Portable.exe" program to open it. You can then use my package normally, for example, using the "ExtSpecR::ExtSpecR_app()" command.

The advantage of using the R-portable compressed package is that it does not require installation, so you can use it directly without affecting the existing system environment. At the same time, all the required packages have been pre-installed, which can save you time and energy and make it easier for you to use my package.

I hope you can use my R package smoothly! If you still encounter any problems, please do not hesitate to contact me and I will do my best to provide assistance.
 
## brief introduction

let me give a quick rundown of our ExtSpecR_app() website.![Screen](/images/first.jpg) We've got three main sections: Introduction, VIs and Examples, and Tree Phenotyping.

In the Introduction section, we've included a handy flowchart that will guide you through the steps of using our app. It's a great way to get a sense of what you can expect and how to get started.

Moving on to VIs and Examples, we have two sub menus: VIs generation and example data. The first one is where you can upload your own five-band TIF raster images from multispectral sensors and calculate different vegetation indices, all while visualizing them on the image display page. Just select the index you want to see, and go to the image display page, You'll get a false-color image and corresponding vegetation index map. Be patient though, it might take some time to process the data due to the size of the TIF files. Once you're done, you can download the resulting images with the click of a button.


The workflow for the VIs calculation is shown in the flowing graph:
![Screenshot](/images/VIs.png)

The page is showing like this:
![Screenshot](/images/viss.png)

In the example data section, we've included some data from my own research, extracted with my ExtSpecR package, which includes tree heights, crown area, and crown layer vegetation indices from research site. You can explore each tree's data by viewing the images and checking out the data info page.

 
Next is the most important part of this app, tree phenotyping. It includes two sub-menus: Data process and Download. The main focus is on the Data process menu, which includes four pages: Upload and Draw ROI, Segmentation, LAS plot information, and Extraction and Visualization. ![Screenshot](/images/treephno.png)

On the Upload and Draw ROI page, users need to upload point cloud data and corresponding raster images data. After uploading, click on the Draw ROI polygon button, and a separate R graphic page will appear. Users can select the area of interest by clicking with the left mouse button and finishing by clicking the right mouse button. The uploaded point cloud data, ROI-selected point cloud data, and image spectral information will be displayed on the Point cloud information and Spectral information sub-page.

After selecting the ROI, users can optimize the Segmentation process by adjusting the hmin and ws information on the Segmentation page![Screenshot](/images/seg.png). After selecting the settings, click on the Start segmentation button, and a small window will appear in the lower right corner of the page, indicating that the calculation is in progress. When the window disappears, two images will appear, a false-color image and the polygon of each tree extracted from the selected ROI area. Users can adjust the TreeID number to locate different tree positions.

After the Segmentation process is complete, users can use the LAS plot information page![Screenshot](/images/las.png) to plot the point cloud information, including the original point cloud data, ROI point cloud data, and segmented trees. The plotted single tree is the one with the TreeID number selected during Segmentation. Users can adjust the point size slider to plot information about different-sized point clouds.

The final step is Extraction and Visualization![Screenshot](/images/exdt.png). After the Segmentation process is complete, users can click on the Extract data button on the Extraction and Visualization page to extract image spectral data on a tree-by-tree basis. Once extracted, the page displays the different data visualizations of each tree's crown, including Z tree height, area crown width, and image spectral information. Users can adjust the Height limits slider to filter the crown spectral information under different height conditions.

Finally, after completing all the steps, users can click on the Download sub-menu to download the final data![Screenshot](/images/down.png), which is available as ".rds" data that can be read in R language. Users can also download the polygon data of each tree extracted during Segmentation. The right-side page shows the data format in an intuitive way.


## data

please download the example data from the bellowed link:[example data](https://www.dropbox.com/sh/dncqmm0eh7ek7sw/AADgg3bgyHGz5HWa-I9wLQxra?dl=0)
or [example2](https://ln5.sync.com/dl/d6899c6f0/3g32725x-b85yuvm3-ba68kfre-jewun6fk)

