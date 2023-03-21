

#' raster_alin aline single raster spectral images together as list
#'
#' @param url the link that contain all raster images
#'
#' @return a list of raster images, used for multi_rasl function
#' @export
#'
#' @examples #please download the example data from this link:https://ln5.sync.com/dl/d6899c6f0/3g32725x-b85yuvm3-ba68kfre-jewun6fk
#' urll1 <- 'D:/Desktop/examples/'
#' images <- raster_alin(urll1)%>% unlist(recursive = F)%>% do.call(c,.)
#' images

#'
#'
#'
#'
raster_alin  <- function(url) {
  lapply(url, function(urll){
    imag <- list.files(
      path = urll,
      pattern = '.tif',
      all.files = T,
      full.names = T,
      no.. = T
    )
    imag <- list(imag)
    # imag1 <-imag[-c(2,9)]
    lapply(imag, function(z)
      expr <- tryCatch({
        library(raster)
        tictoc::tic()
        p_dsm <- raster(z[grepl('dsm|Dsm', z)][1])
        p_blue <- raster::raster(z[grepl('Blue|blue', z)])
        p_green <- raster::raster(z[grepl('Green|green', z)])
        p_red <- raster::raster(z[grepl('Red|red', z)][1])
        p_redege <- raster::raster(z[grepl('dge|dge', z)])
        p_nir <-raster:: raster(z[grepl('NIR|nir', z)])

        p_dsm <- projectRaster(p_dsm, crs='+proj=longlat +datum=WGS84 +no_defs')
        p_dsm <- terra::resample(p_dsm,p_blue, method = 'ngb')
        # extent(p_dsm) <- extent(p_blue)
        # rc1 <- extend(p_nir,p_dsm )
        # p_dsm <- raster::resample(p_dsm, p_nir,method = 'ngb')

        # rrrr <- gdal_resample(p_dsm, p_blue)
        # tictoc::toc()

        trainx <- list(p_red,p_blue,p_green,p_redege,p_nir,p_dsm )
        names(trainx) <- c('red',"blue", "green",'redege','nir','dsm')
        tictoc::toc()
        return(trainx)

      },
      error = function(e) {
        message('Caught an error!')
        cat("ERROR :", conditionMessage(e), "\n")
        print(e)
      },
      print(paste("processing Loop", z, sep = "_"))))})}
