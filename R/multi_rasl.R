#' Do the extraction from raster images using las point data
#'
#'
#'multi_rasl is a function that used in the shiny app for data extraction,
#'it first use the point cloud data to do tree detection and tree crown polygon, and then
#'use this polygon to apply on the multispectal image to extract the multispecral
#'of each tree and then save as a data list
#'
#'
#' @param las_list a list of point cloud data
#' @param dsf_list a list of raster image data, mainly are tif image that read as raster
#' @param kwsindice the parameter for find tree
#' @param hmin  the parameter  of the minimal tree high to use for tree detection. hight less than this values will be abandoned.
#' @return a data list that contain the extracted spectral with x,y, and,treeID,tree hegiht(Z),and crown area
#' @export
#'
#' @examples #please download the example data from this link:https://ln5.sync.com/dl/d6899c6f0/3g32725x-b85yuvm3-ba68kfre-jewun6fk
#'library(lidR)
#'library(ExtSpecR)
#'library(tidyverse)
#'las <- readLAScatalog('D:/Desktop/examples/cloud.las')
#'las$type <- 'example'
#'las_list <- list(las)
#'names(las_list)<- c('example')
#'image_list <- list(images)
#'names(image_list)<- c('example')
#'exam_data <- multi_rasl(las_list,dsf_list =image_list,kwsindice = 7,hmin = 3 )
#'exam_data1 <- exam_data %>% invoke(cbind,.) %>% drop_na()
#'exam_data1[1:5,]
#'
#'
multi_rasl <- function(las_list, dsf_list, kwsindice, hmin) {
  lapply(las_list, function(ctg) {
    expr <- tryCatch({
      library("lidR")
      library("rgdal")
      library(sfheaders)
      library(stars)

      library(tidyverse)
      library(sf)
      library(data.table)
      tictoc::tic("processing las file")
      tictoc::tic("processing rasterize_terrain")
      opt_output_files(ctg) <-
        paste0(tempdir(),rnorm(1), kwsindice, hmin , "{*}_hd")
      opt_chunk_size(ctg) <- 0
      opt_chunk_buffer(ctg) <- 10
      classified_ctg <- classify_ground(ctg, csf())
      dtm <- rasterize_terrain(classified_ctg, 0.5, tin())
      tictoc::toc()
      tictoc::tic("processing normalize_height")
      ctg_norm <- normalize_height(classified_ctg, dtm)
      opt_select(ctg_norm) <- "xyz"
      opt_filter(ctg_norm) <- "-keep_first"
      tictoc::toc()
      tictoc::tic("processing locate_trees")
      opt_output_files(ctg_norm) <- ''
      ttops <- locate_trees(ctg_norm, lmf(ws = kwsindice , hmin = hmin))
      chm <- rasterize_canopy(ctg_norm, 0.5, p2r(0.15))
      tictoc::toc()
      tictoc::tic("processing segment_trees")
      opt_output_files(ctg_norm) <- paste0(tempdir(),rnorm(1), kwsindice, hmin)
      algo <- dalponte2016(chm, ttops)
      ctg_segmented <- segment_trees(ctg_norm, algo)
      tictoc::toc()

      tictoc::tic("processing crown_metrics")
      opt_output_files(ctg_segmented) <- ''
      crown_polo = crown_metrics(ctg_segmented, func = .stdtreemetrics, geom = "convex")
      tictoc::toc()
      directions_longlat <- st_transform(crown_polo, '+proj=longlat +datum=WGS84 +no_defs')
      tictoc::toc()
      tictoc::tic("processing chm2")
      sf  <- st_as_sf(directions_longlat)
      sf33 <- as.data.frame(crown_polo)
      sf33 <-sf33 %>% dplyr::select(treeID, Z, convhull_area) %>% as.data.frame()
      dtm2 <-terra::project(dtm, '+proj=longlat +datum=WGS84 +no_defs')
      chm2 <- terra::project(chm, '+proj=longlat +datum=WGS84 +no_defs')
      tictoc::toc()
      library(data.table)
      library(terra)
      tictoc::tic("processing resample")
      low <- rast(dsf_list[[1]][[1]])
      chm23 <-  terra::resample(chm2, low, method = 'near')
      tictoc::toc()
      tictoc::tic("processing exact_extract")
      library(exactextractr)
      prec_chm <- exactextractr::exact_extract(chm23, sf, include_xy = T) %>%
        setNames(paste0(sf$treeID,"_",sf$convhull_area,"_")) %>%
        purrr::invoke(rbind, .) %>%
        dplyr::select(1:3) %>%
        as.data.frame()
      names(prec_chm)[1] <- 'Z'

      tictoc::toc()
      tictoc::tic("processing spectra exact_extract")
      fd1 <- dsf_list[names(dsf_list) == ctg$type]

      tesst <-  lapply(fd1, function(x11) {
        x22 <- list(unlist(x11))
        lapply(x22, function(ls11) {
          lapply(ls11, function(ls222) {
            tryCatch({
              library(tictoc)
              tic("for loop start")
              message(paste0(names(ls222)))
              print(ls222)
              sp::plot(ls222)
              sp::plot(sf,
                       add = T,
                       alpha = 0.6,
                       col = rainbow(1))
              library(exactextractr)
              prec_dfs <- exactextractr::exact_extract(ls222, sf, include_xy = T) %>%
                setNames(paste0(sf$treeID,"_",sf$convhull_area,"_")) %>%

                purrr::invoke(rbind, .) %>%
                dplyr::select(1) %>%
                as.data.frame()
              names(prec_dfs)  <- names(ls222)
              print("finished")
              toc()
              return(prec_dfs)
            },
            error = function(e) {
              message('Caught an error!')
              cat("ERROR :", conditionMessage(e), "\n")
              print(e)
            })
          })
        })
      }) %>% unlist(recursive = F)  %>% as.data.frame()
      tictoc::toc()

      tictoc::toc()
      rownames(tesst)

      tesst <- tesst[rownames(prec_chm), ]
      dat_tes <- cbind(tesst, prec_chm)
      dat_tes <-  dat_tes %>% mutate(treeID =  sapply(strsplit(rownames(dat_tes), '[_]'), function(x) {
        y = x[1]}))
      dat_tes <-  dat_tes %>% mutate(area =  sapply(strsplit(rownames(dat_tes), '[_]'), function(x) {
        y = x[2]}))
      return(dat_tes)

    }, error = function(e) {
      message('Caught an error!')
      cat("ERROR :", conditionMessage(e), "\n")
      print(e)
    },
    print(paste(
      "processing Loop", names(las_list), sep = "_"
    )))
  })
}
