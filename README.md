
# SpecEXR
![Screen](/images/22.png)

<!-- badges: start -->
<!-- badges: end -->

The SpecEXR is a shiny app for single tree multi-spectral extraction This is a shiny app specially used to extract spectral information of a single tree, and I also provide a sample to show what the data looks like after extraction.
First of all, you need point cloud data with precise positioning information, which is used to segment each individual plant in a large area of forest land. In addition, you need to have multi-spectral or hyperspectral information of this forest land.


## Installation

You can install the development version of SpecEXR from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Yanjie-Li/SpecEXR")
```

Note:`rtools` and "EBImage"  and other packages were required, `rtools` should install manually from: [rtools](https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html) .  
 
 
# The source package

if there are something wrong and can not install from github, please try to download the source packages and install it from your R or Rstudio from package Archive file: [package](/source-package/SpecEXR_1.0.tar.gz)




## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SpecEXR)
## basic example code
SpecexR_app()

```
The first page is "Segmentation", it is an visualization  example that the crown spectral of every single tree on one  plantation site has been extracted. the page is showing like this:
![Screen](/images/segeme2.gif)

The second page is "VIs graphs generation", it is to calculation the VIs from the images based on the five bands raster images, including red, green, blue, rededge and near-infrared raster tif images. 
The workflow for the VIs calculation is shown in the flowing graph:
![Screenshot](/images/VIs.png)



The page is showing like this:
![Screenshot](/images/figure2.png)

The third page is "Tree identification and spectra extraction", it is is the core function of this app, The workflow for how to use this page can be found in the flowing graph:
![Screenshot](/images/treese.png)

The page is showing like this:
![Screenshot](/images/figr23.png)

## data

please download the example data from the bellowed link:[example data](https://www.dropbox.com/sh/dncqmm0eh7ek7sw/AADgg3bgyHGz5HWa-I9wLQxra?dl=0)
or [example2](https://ln5.sync.com/dl/d6899c6f0/3g32725x-b85yuvm3-ba68kfre-jewun6fk)

