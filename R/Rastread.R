#' raster_read
#'
#' @param url
#'
#' @return  a list of images
#' @export
#'
#' @examples  raster_read('url')
#'
#'

raster_read  <- function(url) {

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
        # p_dsm <- raster(z[[1]])
        p_blue <- raster(z[[1]])
        p_green <- raster(z[[2]])
        p_red <- raster(z[[4]])
        p_redege <- raster(z[[5]])
        p_nir <- raster(z[[3]])
        # extent(p_dsm) <- extent(p_blue)
        # rc1 <- extend(p_nir,p_dsm )
        # p_dsm <- raster::resample(p_dsm, p_nir,method = 'ngb')
        # tictoc::tic()
        # rrrr <- gdal_resample(p_dsm, p_blue)
        # tictoc::toc()

        trainx <- list(p_red,p_blue,p_green,p_redege,p_nir )
        # names(trainx) <- c('red',"blue", "green",'redege','nir','dsm')
        return(trainx)

      },
      error = function(e) {
        message('Caught an error!')
        cat("ERROR :", conditionMessage(e), "\n")
        print(e)
      },
      print(paste("processing Loop", z, sep = "_"))))})}

###################
