#' SpecexR_app
#'this is an shiny app mainly for multispectral extraction using LAS point cloud data
#'
#'
#' @param ...
#'
#' @return final data with x,y location information
#' @export
#'
#' @examples  none
#'


SpecexR_app <- function() {
library(shiny)
library(shinydashboard)
library(htmltools)
library(shinyjs)
library(shinyWidgets)
library(shinythemes)
  # packages <- c("shinythemes",'shinyjs', 'RCSF','DT',"shinydashboard",'stars',
  #               'sfheaders','sf','exactextractr', 'lidR', "shiny",'tidyverse',
  #               'RStoolbox','viridis', 'rgdal','tictoc',
  #               'raster','rdrop2','tools','rasterVis','data.table',
  #               "librarian","shinydashboard","tictoc",'BiocManager','quickPlot','pacman')
  # new.packages <- packages[!(packages %in% utils::installed.packages()[,"Package"])]
  # if(length(new.packages)) utils::install.packages(new.packages,repos = "https://cloud.r-project.org")


  # tictoc::tic()
  # # Packages loading
  # if (!require("EBImage", quietly = TRUE))
  #   BiocManager::install('EBImage')
  # if (!require("Biobase", quietly = TRUE))
  #   BiocManager::install("Biobase")
  # librarian::shelf(c(packages,'Biobase'),quiet = T)
  # tictoc::toc()

  print('data read time')
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
      lapply(imag, function(z)
        expr <- tryCatch({
          library(raster)
          tictoc::tic()
          p_dsm <-  raster::raster(z[grepl('dsm|Dsm', z)][1])
          p_blue <- raster::raster(z[grepl('Blue|blue', z)])
          p_green <- raster::raster(z[grepl('Green|green', z)])
          p_red <- raster::raster(z[grepl('Red|red', z)][1])
          p_redege <- raster::raster(z[grepl('dge|dge', z)])
          p_nir <-raster:: raster(z[grepl('NIR|nir', z)])
          p_dsm <- projectRaster(p_dsm, crs='+proj=longlat +datum=WGS84 +no_defs')
          p_dsm <- terra::resample(p_dsm,p_blue, method = 'ngb')
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

  ###################
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
          paste0(tempdir(), kwsindice, hmin , "{*}_hd")
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
        opt_output_files(ctg_norm) <- paste0(tempdir(), kwsindice, hmin)
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
          invoke(rbind, .) %>%
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

                  invoke(rbind, .) %>%
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



  debug_msg <- function(...) {
    is_local <- Sys.getenv('SHINY_PORT') == ""
    in_shiny <- !is.null(shiny::getDefaultReactiveDomain())
    txt <- toString(list(...))
    if (is_local) message(txt)
    if (in_shiny) shinyjs::runjs(sprintf("console.debug(\"%s\")", txt))
  }

  data1 <- readr::read_rds(system.file("data.rds", package = "ExtSpecR"))
  image_path <- system.file('images', "treese.png", package = "ExtSpecR")
  # data1 <- readr::read_rds('inst/data.rds' )
  ras_im_alin <- function(monthi,fami){
    data1 <- readr::read_rds(system.file("data.rds", package = "ExtSpecR"))
    # data1 <- readr::read_rds('inst/data.rds' )
    expr <- tryCatch({
      library(tidyverse)
      library(raster)
      library(EBImage)
      library(tools)
      nir <- filter(data1, month == monthi )
      nir2 <- filter(nir, Fam == fami )
      nir3 <- nir2[,-c(1:4)]
      chl <- nir2[,4]
      library(tidyverse)
      matou_vis <- nir3 %>% dplyr:: mutate(
        ndvi=  ((result_NIR - result_Red) / (result_NIR + result_Red)),
        osavi = ((result_NIR-result_Red)*(1+0.16)) / (result_NIR + result_Red + 0.16),
        gndvi = (result_NIR-result_Green)/(result_NIR+result_Green),
        savi = ((result_NIR - result_Red)*(1+0.5))/((result_NIR + result_Red+0.5)),
        msavi = (2*result_NIR+1-sqrt((2*result_NIR+1)^2-8*(result_NIR-result_Red)))/2,
        gci = result_NIR/result_Green-1,
        RECI = result_NIR/result_RedEdge-1,
        LCI = (result_NIR-result_RedEdge)/(result_NIR+result_Red),
        GRVI =(result_Green-result_Red)/(result_Green+result_Red),
        MGRVI =(result_Green^2-result_Red^2)/(result_Green^2+result_Red^2 ),
        RGBVI =(result_Green^2-result_Red*result_Blue)/(result_Green^2+result_Red*result_Blue),
        NDRE= (result_NIR-result_RedEdge)/(result_NIR+result_RedEdge),
        MACI= result_NIR/result_Green,
        ARI= result_Green/result_NIR,
        MARI=(result_Green^(-1)-result_RedEdge^(-1))/result_NIR

      )

      tryCatch({
        library(raster)
        nir2 <- sapply(matou_vis[,-c(1:2,8:9)], function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))

        matou_vis2 <- cbind.data.frame(matou_vis[, c(1:2)], nir2)
        spg <- matou_vis2
        coordinates(spg) <- ~ x + y
        gridded(spg) <- TRUE
        rasterDF <- stack(spg)
        library(RStoolbox)
        library(rasterVis)
        library(viridis)
        df <-  ggRGB(rasterDF, r=1, g=5, b=3, stretch = 'lin')+ggtitle(paste0(monthi,"_", fami ))
        print(df)
        imgg <-raster:: as.array(rasterDF )
        ds3 <- EBImage::as.Image(imgg)
        ds3 <- ds3[,,-c(1:2)]
        ds4 <- resize(ds3,30,30)
        f_name  <- list(ds4)
        names(f_name) <- paste0(monthi,'_',fami)
        f_namechl  <- list(chl)
        names(f_namechl) <- paste0(monthi,'_chl_',fami)
        message( paste0(monthi,'_',fami))
        return(list(f_namechl,f_name))
      },error = function(e) {NULL})

    },error = function(e) {
      message('Caught an error!')
      cat("ERROR :", conditionMessage(e), "\n")
      print(e)})

  }



  month <- (unique(data1$month))
  fam <- (unique(data1$Fam))

  library(DT)

 css<- "body {
  font-size: 14px;
  font-family: 'Open Sans', sans-serif;
}

.nav > li > a {
  color:  white !important;
  background-color: #ff8c00 !important;
  font-size: 16px;
}
.nav-tabs > li > a {
  color:  white !important;
  background-color: #ff8c00 !important;
  font-size: 16px;
}



.nav-tabs > li > a:active
.nav-tabs>li>a:hover,
.nav-tabs>li>a:focus
{
   background-color: #ffb3e6 !important;
}

.nav > li > a:active,
.nav > li > a:hover,
.nav > li > a:focus {
   background-color: #ffb3e6 !important;
}


.title {
  color: white;
  font-size: 24px !important;
  background-color: #ff8c00 !important;
}

.nav-tabs > li > a {
  color: #FFFFFF;
  background-color: #1f7a8c;
}

.main-header {
  background-color: #1f7a8c !important;
}

.main-sidebar {
  background-color: #334f70 !important;
}

.content-wrapper,
.right-side {
  background-color: #f5f5f5 !important;
}

.content-wrapper {
  padding-top: 15px !important;
}

.sidebar-menu > li > a {
  color: #c7d0d9 !important;
  font-weight: bold !important;
  font-size: 16px !important;
}

.sidebar-menu > li.active > a {
  background-color: #1c4f6d !important;
  color: #c2dd34 !important;
  font-size: 16px !important;
}

.nav-tabs > li > a {
  background-color: #1c4f6d;
  color: #fff;
  font-size: 20px;
}

.treeview-menu > li > a {
  background-color: #1c4f6d;
  color: #c7d0d9 !important;
  font-weight: bold !important;
  font-size: 13px !important;
}

.navbar-default .navbar-brand {
  color: #FFFFFF;
  background-color: #1f7a8c !important;
}

.treeview-menu > li.active > a {
  background-color: #1f7a8c !important;
}

.tab-content {
  background-color: #F7F7F7;
}

.small-box {
  background-color: #fff !important;
  border-radius: 2px !important;
  border: 1px solid #d2d6de !important;
}

.small-box .icon {
  border-radius: 50% !important;
  padding: 10px !important;
  font-size: 24px !important;
  color: #1f7a8c !important;
}

.small-box h3 {
  font-size: 28px !important;
  font-weight: bold !important;
  margin: 0 !important;
}

.small-box p {
  font-size: 14px !important;
  margin: 0 !important;
}
"




# Define the UI
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "Forestry Phenomics",
    titleWidth = "400px"
  ),
  dashboardSidebar(
    sidebarMenu(
      id = "sidebar",
      menuItem(
        "Introduction",
        tabName = "dashboard",
        icon = icon("dashboard")
      ),
      #第一个
      menuItem(
        "VIs and examples",
        icon = icon("database"),
        startExpanded = F,
        menuSubItem(
          "VIs generation",
          tabName = "upload",
          icon = icon("database")
        ),
        menuSubItem(
          "Example data",
          tabName = "edit",
          icon = icon("database")
        )

        ),



   ###第二个
      menuItem(
        "Tree phenotyping",
        icon = icon("database"),
        startExpanded = F,
        # menuSubItem(
        #   "Upload",
        #   tabName = "upload"
        # ),

        menuSubItem(
          "Data input",
          tabName = "input"
        ),
        menuSubItem(
          "Dowload",
          tabName = "calculation"
        )
      )

      # menuItem(
      #   "Settings",
      #   icon = icon("gear"),
      #   startExpanded = FALSE,
      #   menuSubItem(
      #     "General",
      #     tabName = "general"
      #   ),
      #   menuSubItem(
      #     "User",
      #     tabName = "user"
      #   ),
      #   menuSubItem(
      #     "About",
      #     tabName = "about"
      #   )
      # ),

    )
  ),
  dashboardBody(
    includeCSS(system.file("style.css", package = "ExtSpecR")),
    # includeCSS("www/style.css"),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    # tags$head(tags$style(HTML(css))),

    useShinyjs(),
    tabItems(
      tabItem(
        tabName = "dashboard",
        fluidRow(
          box(
            title = "flowchat",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            # height = "500px",
            width = 12,
            fluidPage(

              tags$img(src= system.file("images/treese.png", package = "ExtSpecR") ,width = "100%", height = "auto")
              # 在src中，你需要提供你想要引用的图片的url或者本地路径
              # image_path system.file("images", "treese.png", package = "ExtSpecR")


            )
          )
        )
      ),

      tabItem(
        tabName = "edit",
        navbarPage(
          "Example of tree Segmentation",
          tabPanel(
            "Graphs",
            fluidRow(
              box(
                title = "Plot",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = 12,
                height = "500px",
                sidebarLayout(
                  sidebarPanel(width = 4,
                               selectInput("status", label = h3("fam data"),
                                           choices = fam,
                                           selected = fam[1], multiple = F),
                               selectInput("status2", label =h3("Months data"),
                                           choices = month,
                                           selected = month[1], multiple = F),
                               fluidRow(
                                 numericInput("width_png","Width of PNG", value = 1600) ,
                                 numericInput("height_png","Height of PNG", value = 1200 ),
                                 numericInput("resolution_PNG","Resolution of PNG", value = 144 ),
                                 style = "margin-top: 25px;",
                                 downloadButton('downloadPlotPNG','Download Layer plot PNG'),
                                 downloadButton('downloadPlotPNG2','Download singletree PNG')
                               ),
                               plotOutput('distPlot',width = "100%", height = "800px")



                  ),

                  mainPanel(
                    plotOutput('predictPlot',width = "100%", height = "1000px")

                  )
                )


              )
            )
          ),
          tabPanel(
            "Data Info",
            fluidRow(
              box(
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = 12,
                height = "500px",
                box(title = "Image info", width =6,solidHeader = TRUE,
                     collapsible = TRUE,status = "info", height = "200px", verbatimTextOutput("summary")),
                box(title = "Data info",width =6, solidHeader = TRUE,
                    collapsible = TRUE,status = "info", height = "200px",verbatimTextOutput("summary2"))


              )
            )
          )
        )
      ),



      tabItem(
        tabName = "upload",
        navbarPage(
          "Upload Images",
          tabPanel(
            "Data info",
            fluidRow(
              box(
                title = "please upload related TIF images and select the VIs",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = 12,
                height = "500px",
                fluidRow(
                  column(
                    width = 12,

                    br(),
                    fluidRow(
                      column(
                        width = 6,
                        fileInput("red1", "Upload Red Images", multiple = TRUE, accept = c("tif", ".tif")),
                        fileInput("redege1", "Upload Rededge Images", multiple = TRUE, accept = c("tif", ".tif"))

                      ),
                      column(
                        width = 6,
                        fileInput("gree1", "Upload Green Images", multiple = TRUE, accept = c("tif", ".tif")),
                      ),
                      column(
                        width = 6,
                        fileInput("blue1", "Upload Blue Images", multiple = TRUE, accept = c("tif", ".tif"))
                      ),
                      column(
                        width =12,
                        fileInput("NIR1", "Upload NIR Images", multiple = TRUE, accept = c("tif", ".tif"))
                      )
                    ),
                    pickerInput("testh2o",
                                "VIs:",
                                choices = c(
                                  "Red",
                                  "Green",
                                  "Blue",
                                  "Rededage",
                                  "NIR",
                                  "ndvi",
                                  "osavi",
                                  "gndvi",
                                  "savi",
                                  "msavi",
                                  "gci",
                                  "RECI",
                                  "LCI",
                                  "GRVI",
                                  "MGRVI",
                                  "NDRE",
                                  "MACI",
                                  "ARI",
                                  "MARI"
                                ),
                                selected = "ndvi",
                                multiple = FALSE,
                                options = list(`actions-box` = TRUE),
                                # class = "my-picker-input"
                    )


                  ))
              )
            )
          ),
          tabPanel(
            "Image display",
            column(
              6,
              downloadButton(
                'downloadPlotPNG22',
                'Download RGB TIF',
                # class = "my-download-button"
              )
            ),
            column(
              6,
              downloadButton(
                'downloadPlotPNG11',
                'Download single layer TIF',
                # class = "my-download-button"
              )
            ),

            fluidRow(
              box(
                title = "Images",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = 12,
                height = "500px",

                box(title = "False color image",width =6, solidHeader = TRUE,
                    collapsible = TRUE,status = "info", height = "200px",

                    plotOutput("plotgraph1", width = "100%")),
                box(title = "VI image",width =6, solidHeader = TRUE,
                    collapsible = TRUE,status = "info", height = "200px",

                    plotOutput("plotgraph2", width = "100%"))

              )


            )
          )
        )
      ),

      tabItem(
        tabName = "input",
        navbarPage(
          "View Data",
          tabPanel(

            "Upload and Draw ROI",
            fluidRow(
              box(
                title = "Upload data and start draw ROI",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = 6,
                height = "500px",



                fluidRow(
                  column(6,
                         fluidRow(
                           column(12, "Point cloud File"),
                           column(12,
                                  tags$div(style = "background-color: transparent;",
                                    fileInput("file1", "", multiple = FALSE, accept = c("las", ".las", ".laz", ".ply"))
                                    # tags$p("Choose a point cloud file")
                                  )
                           )
                         )
                  ),
                  column(6,
                         fluidRow(
                           column(12, "Raster File"),
                           column(12,
                                  tags$div(style = "background-color: transparent;",
                                    fileInput("file2", "", multiple = TRUE, accept = c("tif", ".tif"))
                                    # tags$p("Choose a raster file")
                                  )
                           )
                         )
                  )
                ),
                fluidRow(
                  column(12,
                         tags$div(  style = "text-align: center;", # 将内容居中,
                           actionButton("drawpoly", "Draw ROI polygon", class = "primary"),
                           tags$p("Click the button to draw an ROI region"),

                           tags$head(tags$style("#drawpoly {
  color: #fff;
  font-size: 50px;
  font-weight: bold;
  text-shadow: 1px 1px #888888;
  background-color: #0B658A;
  border: none;
  border-radius: 5px;
  box-shadow: 0px 3px 2px rgba(0, 0, 0, 0.2);
}

#drawpoly:hover {
  background-color: #135E8C;
  box-shadow: 0px 6px 5px rgba(0, 0, 0, 0.2);
  transform: translateY(-2px);
}

#drawpoly:active {
  box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.2);
  transform: translateY(1px);
}"
                           ))
                         )
                  )
                )


              ),
              box(
                title = "information",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = 6,
                height = "500px",
                h1("Briefly introduction"),
                p("Please use cloud point",
                  em("(.las)"),
                  "and raster images",
                  em("(.tif)"),
                  "to start the extraction"),
                tags$hr(),
                fluidRow(
                  tabsetPanel(type = "tabs",
                              tabPanel("Point cloud information",
                                       verbatimTextOutput("summar"),
                                       verbatimTextOutput("sudra")
                              ),
                              tabPanel("Spectral information",
                                       fluidRow(
                                           verbatimTextOutput("summar2")

                                       )
                              )
                  )
                )


              )
            )
          ),

          tabPanel(
            "Segmentation",



            box(title = "Choice parameter and click the start segementation",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = 12,
                height = "500px",
             fluidRow(
              column(12, br()),
              column(4,
                     tags$div(style = "background-color: transparent;", selectInput(
                       inputId = "wscontro",
                       label = 'ws:',
                       choices = seq(0, 10, by = 0.1),
                       selected = 6
                     ))
              ),
              column(4,
                     tags$div(style = "background-color: transparent;",  selectInput(
                       inputId = "hmincor",
                       label = 'hmin:',
                       choices = seq(1, 10, by = 0.1),
                       selected = 1.5
                     ))
              ),
              column(4,
                     tags$div(style = "background-color: transparent;",  selectInput(
                       inputId = "select2",
                       label = "TreeID number:",
                       choices = c(1:5000),
                       selected = 2
                     ))
              ),
              column(12,
                     tags$div(
                       style = "text-align: center;", # 将内容居中
                       actionButton("statdata", "Start segmentation", class = "primary"),
                       tags$p("Click start to do individual tree detection and segmentation")
                     )
              ),

              tags$head(tags$style("#statdata {
  color: #fff;
  font-size: 50px;
  font-weight: bold;
  text-shadow: 1px 1px #888888;
  background-color: #0B658A;
  border: none;
  border-radius: 5px;
  box-shadow: 0px 3px 2px rgba(0, 0, 0, 0.2);
}

#statdata:hover {
  background-color: #135E8C;
  box-shadow: 0px 6px 5px rgba(0, 0, 0, 0.2);
  transform: translateY(-2px);
}

#statdata:active {
  box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.2);
  transform: translateY(1px);
}"
              ))


            ),



            fluidRow(
              column(4,  tags$div(
                style = "text-align: center;", # 将内容居中
                downloadButton('downloadrgball',
                               'Downloadrgb',class = "butt1"))  ),

              tags$head(tags$style("#downloadrgball {
  color: #fff;
  font-size: 15px;
  font-weight: normal;
  text-shadow: 1px 1px #888888;
  background-color: #0B658A;
  border: none;
  border-radius: 5px;
  box-shadow: 0px 3px 2px rgba(0, 0, 0, 0.2);
}

#downloadrgball:hover {
  background-color: #135E8C;
  box-shadow: 0px 6px 5px rgba(0, 0, 0, 0.2);
  transform: translateY(-2px);
}

#downloadrgball:active {
  box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.2);
  transform: translateY(1px);
}"
              )),

              column(4,  tags$div(
                style = "text-align: center;", # 将内容居中
                downloadButton('downloadsfall','Downloadshp',class = "butt1"))  ),

              tags$head(tags$style("#downloadsfall {
  color: #fff;
  font-size: 15px;
  font-weight: normal;
  text-shadow: 1px 1px #888888;
  background-color: #0B658A;
  border: none;
  border-radius: 5px;
  box-shadow: 0px 3px 2px rgba(0, 0, 0, 0.2);
}

#downloadsfall:hover {
  background-color: #135E8C;
  box-shadow: 0px 6px 5px rgba(0, 0, 0, 0.2);
  transform: translateY(-2px);
}

#downloadsfall:active {
  box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.2);
  transform: translateY(1px);
}"
              ))

            ), fluidRow(
              verbatimTextOutput("warnpredict") ,
              tags$hr(),

              verbatimTextOutput("warnpredict2"),
              fluidRow(  splitLayout(
                style = "border: 1px solid silver:",

                plotOutput("predictPlo" , width = "100%", height = "800px" ),

                plotOutput("predictPlot5", width = "100%", height = "800px")
              ))



            ) )

          ),
          tabPanel(
            "LAS polt information",
            fluidRow(
              box(
                title = "LAS polt information",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = 12,
                height = "500px",


                fluidRow(
                  column(3,  tags$div(
                    style = "text-align: center;", # 将内容居中
                    actionButton("dodo", "Plot las point cloud"))  ),

                  tags$head(tags$style("#dodo {
  color: #fff;
  font-size: 20px;
  font-weight: bold;
  text-shadow: 1px 1px #888888;
  background-color: #0B658A;
  border: none;
  border-radius: 5px;
  box-shadow: 0px 3px 2px rgba(0, 0, 0, 0.2);
}

#dodo:hover {
  background-color: #135E8C;
  box-shadow: 0px 6px 5px rgba(0, 0, 0, 0.2);
  transform: translateY(-2px);
}

#dodo:active {
  box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.2);
  transform: translateY(1px);
}"
                  )),




                  column(3,  tags$div(
                    style = "text-align: center;", # 将内容居中
                    actionButton("drawpolyroi", "Plot ROI point cloud"))  ),

                  tags$head(tags$style("#drawpolyroi {
  color: #fff;
  font-size: 20px;
  font-weight: bold;
  text-shadow: 1px 1px #888888;
  background-color: #0B658A;
  border: none;
  border-radius: 5px;
  box-shadow: 0px 3px 2px rgba(0, 0, 0, 0.2);
}

#drawpolyroi:hover {
  background-color: #135E8C;
  box-shadow: 0px 6px 5px rgba(0, 0, 0, 0.2);
  transform: translateY(-2px);
}

#drawpolyroi:active {
  box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.2);
  transform: translateY(1px);
}"
                  )),

                  column(3,  tags$div(
                    style = "text-align: center;", # 将内容居中
                    actionButton("dodo1", "Plot single tree point cloud"))  ),



                  tags$head(tags$style("#dodo1 {
  color: #fff;
  font-size: 20px;
  font-weight: bold;
  text-shadow: 1px 1px #888888;
  background-color: #0B658A;
  border: none;
  border-radius: 5px;
  box-shadow: 0px 3px 2px rgba(0, 0, 0, 0.2);
}

#dodo1:hover {
  background-color: #135E8C;
  box-shadow: 0px 6px 5px rgba(0, 0, 0, 0.2);
  transform: translateY(-2px);
}

#dodo1:active {
  box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.2);
  transform: translateY(1px);
}"
                  )),


                  column(12, sliderInput(inputId = 'poinsize',
                                         label = 'point size:',
                                         value = 2,
                                         step =1,
                                         min =0 ,
                                         max =10)),

                  fluidRow(

                    column(3,  tags$div(
                      style = "text-align: center;", # 将内容居中
                      downloadButton("tif_data", label = "Download spectral data"))  ),

                    tags$head(tags$style("#tif_data {
  color: #fff;
  font-size: 15px;
  font-weight: bold;
  text-shadow: 1px 1px #888888;
  background-color: #0B658A;
  border: none;
  border-radius: 5px;
  box-shadow: 0px 3px 2px rgba(0, 0, 0, 0.2);
}

#tif_data:hover {
  background-color: #135E8C;
  box-shadow: 0px 6px 5px rgba(0, 0, 0, 0.2);
  transform: translateY(-2px);
}

#tif_data:active {
  box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.2);
  transform: translateY(1px);
}"
                    )),



                    column(3,  tags$div(
                      style = "text-align: center;", # 将内容居中
                      downloadButton("cleaned_data", label = "Download ROI PCD (las)"))  ),

                    tags$head(tags$style("#cleaned_data {
  color: #fff;
  font-size: 15px;
  font-weight: bold;
  text-shadow: 1px 1px #888888;
  background-color: #0B658A;
  border: none;
  border-radius: 5px;
  box-shadow: 0px 3px 2px rgba(0, 0, 0, 0.2);
}

#cleaned_data:hover {
  background-color: #135E8C;
  box-shadow: 0px 6px 5px rgba(0, 0, 0, 0.2);
  transform: translateY(-2px);
}

#cleaned_data:active {
  box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.2);
  transform: translateY(1px);
}"
                    )),


                    column(3,  tags$div(
                      style = "text-align: center;", # 将内容居中
                      downloadButton("cleaned_data2", label = "Download ROI PCD (laz)"))  ),

                    tags$head(tags$style("#cleaned_data2 {
  color: #fff;
  font-size: 15px;
  font-weight: bold;
  text-shadow: 1px 1px #888888;
  background-color: #0B658A;
  border: none;
  border-radius: 5px;
  box-shadow: 0px 3px 2px rgba(0, 0, 0, 0.2);
}

#cleaned_data2:hover {
  background-color: #135E8C;
  box-shadow: 0px 6px 5px rgba(0, 0, 0, 0.2);
  transform: translateY(-2px);
}

#cleaned_data2:active {
  box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.2);
  transform: translateY(1px);
}"
                    ))
                  )
                  ,

                  fluidRow(
                    column(4, numericInput("width_png3","Width of PNG", value = 1600)) ,
                    column(4, numericInput("height_png3","Height of PNG", value = 1200 )),
                    column(4, numericInput("resolution_PNG3","Resolution of PNG", value = 144 )),
                    column(4, numericInput("width_pdf3","Width of pdf", value = 16)) ,
                    column(4, numericInput("height_pdf3","Height of pdf", value = 12 )),

                    style = "margin-top: 25px;")
                ),
                fluidRow(  splitLayout(
                  style = "border: 1px solid silver:",

                  plotOutput("contents" , width = "100%", height = "800px"),

                  plotOutput("plotfinalresult", width = "80%", height = "800px")
                )) ,

                tableOutput("contents33"),
                tableOutput("contents22")





              )
            )
          ),

          tabPanel(
            "Extraction and visualization",
            fluidRow(
              box(
                title = "Extraction and visualization",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = 12,
                height = "500px",

                column(6,
                       tags$div(
                         style = "text-align: center;", # 将内容居中
                         actionButton("extractBtn", "Extract data", class = "primary"),
                         tags$p("Extract individual tree Spectral")
                       )
                ),

                tags$head(tags$style("#extractBtn {
      color: #fff;
      font-size: 50px;
      font-weight: bold;
      text-shadow: 1px 1px #888888;
      background-color: #0B658A;
      border: none;
      border-radius: 5px;
      box-shadow: 0px 3px 2px rgba(0, 0, 0, 0.2);
    }

    #extractBtn:hover {
      background-color: #135E8C;
      box-shadow: 0px 6px 5px rgba(0, 0, 0, 0.2);
      transform: translateY(-2px);
    }

    #extractBtn:active {
      box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.2);
      transform: translateY(1px);
    }"
                )),

                fluidRow(
                  column(6, sliderInput("heightdata",
                                        "Height limits:",
                                        min = 0,
                                        max = 10,
                                        step = 0.1,
                                        value = 0))

                ),

                fluidRow(column(12, splitLayout(
                  style = "border: 1px solid silver;",
                  plotOutput("predictPlot3" , width = "100%", height = "800px"),
                  plotOutput("predictPlot4", width = "100%", height = "800px")
                )))
              )
            )

          )





        )
      ),

      tabItem(
        tabName = "calculation",
        navbarPage(
          "Download",
          tabPanel(
            "Final data",
            fluidRow(
              box(
                # title = "about Data - Tab 1",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = 6,
                height = "500px",

                fluidRow(
                  tags$div(
                    style = "display: flex; justify-content: center;",
                    column(6, downloadButton("tife_data", label = "Download final extracted data")),
                    column(6, downloadButton("sf_data", label = "Download final shapefile data"))
                  ),

                  tags$head(tags$style(
                    "#tife_data, #sf_data {
      color: #fff;
      font-size: 20px;
      font-weight: bold;
      text-shadow: 1px 1px #888888;
      background-color: #0B658A;
      border: none;
      border-radius: 5px;
      box-shadow: 0px 3px 2px rgba(0, 0, 0, 0.2);
    }

    #tife_data:hover, #sf_data:hover {
      background-color: #135E8C;
      box-shadow: 0px 6px 5px rgba(0, 0, 0, 0.2);
      transform: translateY(-2px);
    }

    #tife_data:active, #sf_data:active {
      box-shadow: 0px 2px 2px rgba(0, 0, 0, 0.2);
      transform: translateY(1px);
    }"
                  )),

                  tags$hr()
                )


              ),
              box( status = "primary",
                   solidHeader = TRUE,
                   collapsible = TRUE,
                   width = 6,
                   height = "500px",

              mainPanel(verbatimTextOutput("su2"))
              )
            )
          )
        )
      )

    )
  )
)
# 服务器部分


library(tools)
options(shiny.maxRequestSize=5000*1024^2)
server <- function(input, output) {
  trend_data <- reactive({
    sele <-input$status
  })
  trend_data2 <- reactive({
    selec <-input$status2
  })

  output$summary <- renderPrint({
    chl_tree <- sigletree2()
    gfdgh1 <- chl_tree %>%unlist(recursive = F)
    gfdgh1[sapply(gfdgh1, is.null)] <- NULL
    library(data.table)
    fffghj2 <- (gfdgh1[!names(gfdgh1) %like% 'chl'])
    fffghj2[sapply(fffghj2, is.null)] <- NULL
    fffghj2[[1]][is.na(fffghj2[[1]])] <- 0
    print(fffghj2[[1]])

  })



    if (!interactive()) {
      session$onSessionEnded(function() {
        stopApp()
        q("no")
      })
    }






  sigletree <- reactive({

    library(tidyverse)
    library(raster)
    library(EBImage)
    library(tools)

    nir <- filter(data1, month == trend_data2() )
    nir2 <- filter(nir, Fam == trend_data() )
    nir3 <- nir2[,-c(1:4)]
    chl <- nir2[,4]
    library(tidyverse)
    matou_vis <- nir3 %>% dplyr:: mutate(
      ndvi=  ((result_NIR - result_Red) / (result_NIR + result_Red)),
      osavi = ((result_NIR-result_Red)*(1+0.16)) / (result_NIR + result_Red + 0.16),
      gndvi = (result_NIR-result_Green)/(result_NIR+result_Green),
      savi = ((result_NIR - result_Red)*(1+0.5))/((result_NIR + result_Red+0.5)),
      msavi = (2*result_NIR+1-sqrt((2*result_NIR+1)^2-8*(result_NIR-result_Red)))/2,
      gci = result_NIR/result_Green-1,
      RECI = result_NIR/result_RedEdge-1,
      LCI = (result_NIR-result_RedEdge)/(result_NIR+result_Red),
      GRVI =(result_Green-result_Red)/(result_Green+result_Red),
      MGRVI =(result_Green^2-result_Red^2)/(result_Green^2+result_Red^2 ),
      RGBVI =(result_Green^2-result_Red*result_Blue)/(result_Green^2+result_Red*result_Blue),
      NDRE= (result_NIR-result_RedEdge)/(result_NIR+result_RedEdge),
      MACI= result_NIR/result_Green,
      ARI= result_Green/result_NIR,
      MARI=(result_Green^(-1)-result_RedEdge^(-1))/result_NIR

    )

    tryCatch({
      library(raster)
      nir2 <- sapply(matou_vis[,-c(1:2,8:9)], function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))

      matou_vis2 <- cbind.data.frame(matou_vis[, c(1:2)], nir2)
      spg <- matou_vis2
      coordinates(spg) <- ~ x + y
      gridded(spg) <- TRUE
      rasterDF <- stack(spg)
      library(RStoolbox)
      library(rasterVis)
      library(viridis)
      df <-  ggRGB(rasterDF, r=1, g=5, b=3, stretch = 'lin')+ggtitle(paste0(trend_data2(),"_",  trend_data() ))
      print(df)
    })

  })
  output$distPlot <- renderPlot({
    chl_tree <- sigletree()
    print(chl_tree)

  })


  sigletree2 <- reactive({

    chl_treew <-  ras_im_alin(fami= trend_data(),monthi=trend_data2())
    chl_treew

  })


  output$summary2  <- renderPrint({
    library(tidyverse)
    library(raster)
    library(EBImage)
    library(tools)
    nir <- filter(data1, month == trend_data2() )
    nir2 <- filter(nir, Fam == trend_data() )
    nir3 <- nir2[,-c(1:4)]
    chl <- nir2[,4]
    library(tidyverse)
    matou_vis <- nir3 %>% dplyr:: mutate(
      ndvi=  ((result_NIR - result_Red) / (result_NIR + result_Red)),
      osavi = ((result_NIR-result_Red)*(1+0.16)) / (result_NIR + result_Red + 0.16),
      gndvi = (result_NIR-result_Green)/(result_NIR+result_Green),
      savi = ((result_NIR - result_Red)*(1+0.5))/((result_NIR + result_Red+0.5)),
      msavi = (2*result_NIR+1-sqrt((2*result_NIR+1)^2-8*(result_NIR-result_Red)))/2,
      gci = result_NIR/result_Green-1,
      RECI = result_NIR/result_RedEdge-1,
      LCI = (result_NIR-result_RedEdge)/(result_NIR+result_Red),
      GRVI =(result_Green-result_Red)/(result_Green+result_Red),
      MGRVI =(result_Green^2-result_Red^2)/(result_Green^2+result_Red^2 ),
      RGBVI =(result_Green^2-result_Red*result_Blue)/(result_Green^2+result_Red*result_Blue),
      NDRE= (result_NIR-result_RedEdge)/(result_NIR+result_RedEdge),
      MACI= result_NIR/result_Green,
      ARI= result_Green/result_NIR,
      MARI=(result_Green^(-1)-result_RedEdge^(-1))/result_NIR

    )

    library(raster)
    nir2s <- sapply(matou_vis[,-c(1:2,8:9)], function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))

    matou_vis2 <- cbind.data.frame(nir2[, c(1:6,12:13)], nir2s)
    matou_vis2
  })

  output$predictPlot <- renderPlot({
    chl_tree <- sigletree2()
    gfdgh1 <- chl_tree %>%unlist(recursive = F)
    gfdgh1[sapply(gfdgh1, is.null)] <- NULL
    library(data.table)
    fffghj2 <- (gfdgh1[!names(gfdgh1) %like% 'chl'])
    fffghj2[sapply(fffghj2, is.null)] <- NULL
    fffghj2[[1]][is.na(fffghj2[[1]])] <- 0
    y <- brick(fffghj2[[1]])
    library(viridis)
    sp::plot(y,col=viridis(20))


  })

  df_products_upload <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return('please upload las cloud data')
    las_12 <- lapply(inFile$datapath,function(m){
      fdd <- lidR::readLAScatalog(m )

    } )

    las_12

  })




  adraw <-  reactive({
    df_point <- df_products_upload()
    #
    las <-readLAS(df_point[[1]])
    chm <- rasterize_canopy(las, res = 0.1, p2r())
    x11()
    plot.new()
    sp:: plot(chm,col=height.colors(150))
    hf <- terra::draw(x="polygon", id=T)
    df <- sf::st_as_sf(hf)
    ctg_sd <- df_point[[1]]
    opt_output_files(ctg_sd) <-paste0(tempdir(),rnorm(1),'_{ID}' )
    x11()
    plot.new()
    subset3 <- lidR::clip_roi(ctg_sd,df)

  })


  adraw3 <- eventReactive(input$drawpoly, {
    withProgress(message = 'External window will open',
                 detail = 'please wait...', value = 0, {
                   sele <- adraw()
                 })
  })

  output$sudra  <- renderPrint({
    adraw5 <- adraw3()
    print(adraw5)
  })

  draw_cloud   <- reactive({
    adraw5 <- adraw3()
    adraw5
  })







  polyroi <- eventReactive(input$drawpolyroi, {
    withProgress(message = 'Ploting',
                 detail = 'May take a while...', value = 0, {
                   sele <- draw_cloud()
                   sele <- readLAS(sele)
                   # lapy <- lapply(sele, function(x){
                     sp::plot(sele, bg = "white",size = input$poinsize  , axis = TRUE, legend = TRUE)
                   # })
                   lapy
                 })
  })




  output$contents33 <- renderPrint({
    withProgress(message = 'Ploting',
                 detail = 'May take a while...',
                 value = 0,
                 {
                   library("lidR")
                   library("rgdal")
                   library(raster)
                   library(tidyverse)
                   print("plot with RGL device")
                   rma <-  polyroi()

                   if (is.null(rma))
                     return(NULL)
                   print(rma)
                 })

  })


  output$downloadPlotPNG <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("ggtree_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {

      png(file, width = input$width_png, height = input$height_png, res = input$resolution_PNG)
      chl_tree <- sigletree2()
      gfdgh1 <- chl_tree %>%unlist(recursive = F)
      gfdgh1[sapply(gfdgh1, is.null)] <- NULL
      library(data.table)
      fffghj2 <- (gfdgh1[!names(gfdgh1) %like% 'chl'])
      fffghj2[sapply(fffghj2, is.null)] <- NULL
      fffghj2[[1]][is.na(fffghj2[[1]])] <- 0
      y <- brick(fffghj2[[1]])
      library(viridis)
      sp::plot(y,col=viridis(20))
      dev.off()
    },

    contentType = "application/png" # MIME type of the image

  )


  output$downloadPlotPNG2 <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("ggtree_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {

      png(file, width = input$width_png, height = input$height_png, res = input$resolution_PNG)
      ras_im_alin(fami= trend_data(),monthi=trend_data2())
      dev.off()
    },

    contentType = "application/png"

  )


  getData <- reactive({
    library(data.table)
    inFile1 <- input$file2
    if (is.null(inFile1)){
      return(print('please upload raster images')) } else{
        withProgress(message = 'Calculation in progress',
                     detail = 'This may take a while...', value = 0, {
                       plr <- lapply(inFile1$datapath, function(m){
                         rs <-terra::rast(m)

                         if( !(terra::crs(rs,  proj=TRUE)) == '+proj=longlat +datum=WGS84 +no_defs' ){
                           message(paste('projectioning'))
                           tictoc::tic(print('projection'))
                           rs2 <- terra::project(rs, '+proj=longlat +datum=WGS84 +no_defs')
                           tictoc::toc()
                           rs2

                         } else{rs}


                       }) %>%    do.call(c,.)

                       plr
                     })
      }

  })


  output$tif_data  <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("spetral_", gsub("/", "-", x), ".tif", sep = "")
    },
    content = function(file) {
      withProgress(message = 'Downloading',
                   detail = 'please wait...', value = 0, {

                     dsf1 <-  getData()
                     dsf1 <-  dsf1
                     res <- terra::writeRaster(dsf1, filename=file,gdal=c("COMPRESS=DEFLATE", "TFW=YES"), overwrite=TRUE, datatype='INT1U')
                   })
    }

  )










  output$cleaned_data <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("ROI_", gsub("/", "-", x), ".las", sep = "")
    },
    content = function(file) {
      withProgress(message = 'Downloading',
                   detail = 'please wait...', value = 0, {
                     sele <- draw_cloud()
                     sele <- readLAS(sele)
                     writeLAS(sele, file, index= TRUE)

                   })
    }


  )

  output$cleaned_data2 <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("ROI_", gsub("/", "-", x), ".laz", sep = "")
    },
    content = function(file) {
      withProgress(message = 'Downloading',
                   detail = 'please wait...', value = 0, {
                     sele <- draw_cloud()
                     sele <- readLAS(sele)
                     writeLAS(sele, file, index= TRUE)

                   })
    }


  )

  output$summar  <- renderPrint({
    las22 <- df_products_upload()
    print(las22)
  })



  data_dsf1 <- reactive({
    select1 <-input$caption
    library("lidR")
    library("rgdal")
    library(raster)
    library(tidyverse)
    dsf1 <- raster_read(select1)
    dsf1 <- dsf1 %>% unlist(recursive = F) %>%  unlist(recursive = F)

  })


  output$summar2  <- renderPrint({
    dsf1 <-  getData()
    if (is.null(dsf1))
      return(NULL)
    print(dsf1)

  })



  mult  <-  reactive({
    las_12 <- draw_cloud()
    las_list <- list(las_12)
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...', value = 0, {
                   lasdata <-   lapply(las_list, function(ctg){
                     expr <- tryCatch({
                       library("lidR")
                       library("rgdal")
                       library(sfheaders)
                       library(stars)
                       library(raster)
                       library(tidyverse)
                       library(sf)
                       library(data.table)

                       tictoc:: tic("processing las file")
                       tictoc:: tic("processing dtm")
                       ctg$overwrite <- TRUE
                       opt_output_files(ctg) <-paste0(tempdir(),rnorm(1), as.numeric( input$wscontro),as.numeric(input$hmincor) ,"/{ORIGINALFILENAME}_{ID}")

                       opt_chunk_size(ctg) <- 0
                       opt_chunk_buffer(ctg) <- 20
                       classified_ctg <- classify_ground(ctg, csf())
                       dtm <- rasterize_terrain(classified_ctg, 1, tin())
                       tictoc:: toc()
                       tictoc:: tic("processing normalize_h")
                       ctg_norm <- normalize_height(classified_ctg, dtm)
                       opt_select(ctg_norm) <- "xyz"
                       opt_filter(ctg_norm) <- "-keep_first"
                       hmean3 <- pixel_metrics(ctg_norm, ~mean(Z), 0.2)
                       tictoc:: toc()
                       hmean3[hmean3 < 0.1] <- NA

                       tictoc:: tic("processing ttops")
                       ctg_norm$overwrite <- TRUE
                       opt_output_files(ctg_norm) <-''
                       ttops <- locate_trees(hmean3,  lmf(ws=as.numeric( input$wscontro) , hmin = as.numeric(input$hmincor)))


                       tictoc:: toc()
                       tictoc:: tic("processing chm")
                       chm <- rasterize_canopy(ctg_norm, 0.2, p2r(0.15))
                       tictoc:: toc()
                       tictoc:: tic("processing segment_trees")
                       opt_output_files(ctg_norm) <-paste0(tempdir(),rnorm(1),as.numeric( input$wscontro),as.numeric(input$hmincor),'{ORIGINALFILENAME}_{XCENTER}_{ID}' )
                       algo <- dalponte2016(chm, ttops)
                       ctg_segmented <- segment_trees(ctg_norm, algo)
                       tictoc:: toc()
                       opt_output_files(ctg_segmented) <- ''
                       crown_polo = crown_metrics(ctg_segmented, func = .stdtreemetrics, geom = "convex")
                       directions_longlat <-  st_transform(crown_polo, '+proj=longlat +datum=WGS84 +no_defs')
                       sf  <- st_as_sf(directions_longlat)
                       sf  <- sf%>% dplyr:: select(treeID,convhull_area)
                       chm2 <- terra::project(chm, '+proj=longlat +datum=WGS84 +no_defs')

                       lidR:::catalog_laxindex(ctg_segmented)
                       tem3 <- list(crown_polo,sf,chm2,ctg_segmented)
                       names(tem3) <- c('crown_polo','sf','chm2','ctg_segmented')
                       return(tem3)
                     },error = function(e) {
                       message('Caught an error!')
                       cat("ERROR :", conditionMessage(e), "\n")
                       print(e)},
                     print(paste("processing Loop", names(las_list), sep = "_")))
                   }) %>%  unlist(recursive = F)

                   lasdata })
  })




  output$sf_data <- downloadHandler(

    filename <- function() {
      "Data_shpExport.zip"
    },
    content = function(file) {
      withProgress(message = "Exporting Data", {
        library(sf)
        library(zip)

        incProgress(0.5)
        tmp.path <- dirname(file)

        name.base <- file.path(tmp.path, "Individual tree shp")
        name.glob <- paste0(name.base, ".*")
        name.shp  <- paste0(name.base, ".shp")
        name.zip  <- paste0(name.base, ".zip")

        if (length(Sys.glob(name.glob)) > 0) file.remove(Sys.glob(name.glob))
        crowte <- adrarr()$crown_polo
        sf::st_write(crowte, dsn = name.shp, ## layer = "shpExport",
                     driver = "ESRI Shapefile", quiet = TRUE)

        zip::zipr(zipfile = name.zip, files = Sys.glob(name.glob))
        req(file.copy(name.zip, file))

        if (length(Sys.glob(name.glob)) > 0) file.remove(Sys.glob(name.glob))

        incProgress(0.5)
      })
    }
  )



  adrarr <- eventReactive(input$statdata, {
    sele <- mult()
    sele
  })

  extractData <- function() {
    dsf1 <- getData()

    if (is.null(dsf1))
      return('please upload raster images')

    withProgress(message = 'Extraction in progress',
                 detail = 'Time consuming...,please wait', value = 0, {
                   crowte  <- adrarr()$crown_polo
                   ctg_segmented <- adrarr()$ctg_segmented

                   def <- st_crs(ctg_segmented)

                   pr <- terra::project(dsf1,def$input,res=0.1)

                   crff  <- 1:length (crowte$treeID)
                   crown_se  <- lapply(crff, function(fdx){
                     expr <- tryCatch({

                       sf2  <-  crowte %>% dplyr:: filter(treeID == fdx)
                       subset3 <- clip_roi(ctg_segmented, sf2)

                       fer <- payload(subset3)  %>% dplyr::select(X,Y,Z,treeID)
                       fer$treeID <- as.factor(fer$treeID)
                       names(fer) <- c('x','y','z','treeID')
                       fer <- as.data.frame(fer)

                       message(paste0('project',fdx))

                       dsta <- terra::extract(pr,fer[,c('x','y')],xy=T ) %>% mutate(treeID=fer$treeID,
                                                                                    Z=fer$z,
                                                                                    area=sf2$convhull_area
                       ) %>%drop_na()
                       return(dsta)

                     },error = function(e){
                       message('Caught an error!')
                       paste(NaN)

                     })
                   }) %>% invoke(rbind,.)  %>% as.data.frame()%>% mutate_if(is.character,as.numeric)

                   return(crown_se)
                 })
  }

  data_ext2 <- eventReactive(input$extractBtn, {
    extractData()
  })


  # data_ext2  <- reactive({
  #   dsf1 <- getData()
  #
  #   if (is.null(dsf1))
  #     return('please upload raster images')
  #
  #   withProgress(message = 'Extraction in progress',
  #                detail = 'Time consuming...,please wait', value = 0, {
  #                  crowte  <- adrarr()$crown_polo
  #                  ctg_segmented <- adrarr()$ctg_segmented
  #
  #                  def <- st_crs(ctg_segmented)
  #
  #                  pr <- terra::project(dsf1,def$input,res=0.1)
  #
  #                  crff  <- 1:length (crowte$treeID)
  #                  crown_se  <- lapply(crff, function(fdx){
  #                    expr <- tryCatch({
  #
  #                      sf2  <-  crowte %>% dplyr:: filter(treeID == fdx)
  #                      subset3 <- clip_roi(ctg_segmented, sf2)
  #
  #                      fer <- payload(subset3)  %>% dplyr::select(X,Y,Z,treeID)
  #                      fer$treeID <- as.factor(fer$treeID)
  #                      names(fer) <- c('x','y','z','treeID')
  #                      fer <- as.data.frame(fer)
  #
  #                      message(paste0('project',fdx))
  #
  #                      dsta <- terra::extract(pr,fer[,c('x','y')],xy=T ) %>% mutate(treeID=fer$treeID,
  #                                                                                   Z=fer$z,
  #                                                                                   area=sf2$convhull_area
  #                      ) %>%drop_na()
  #                      return(dsta)
  #
  #                    },error = function(e){
  #                      message('Caught an error!')
  #                      paste(NaN)
  #
  #                    })
  #                  }) %>% invoke(rbind,.)  %>% as.data.frame()%>% mutate_if(is.character,as.numeric)
  #
  #                  crown_se
  #                })
  #
  # })


  sertree  <-  reactive({
    withProgress(message = 'Ploting',
                 detail = 'May take a while...',
                 value = 0,
                 {
                   crowte  <- adrarr()$crown_polo
                   ctg_segmented <- adrarr()$ctg_segmented
                   sf2  <-
                     crowte %>% dplyr::filter(treeID == as.numeric(input$select2))
                   subset3 <- clip_roi(ctg_segmented, sf2)
                   subset2 <-
                     filter_poi(subset3, Z >= input$heightdata)

                   sp::plot(
                     subset2,
                     bg = "white",
                     size = input$poinsize,
                     axis = TRUE,
                     legend = TRUE
                   )
                 })

  })





  randse <- eventReactive(input$dodo1, {
    sele <- sertree()
  })

  output$contents  <- renderPlot({
    library("lidR")
    library("rgdal")
    library(raster)
    library(tidyverse)
    print("plot with RGL device")
    randese <-  randse()
    if (is.null(randese))
      return(NULL)
    print(randese)

  })
  output$plotfinalresult  <- renderPlot({
    finalda  <- finaldata()
    nir <- filter(finalda, treeID == as.numeric(input$select2) )
    re <- names(nir)
    p1 <-  nir %>% filter(Z> input$heightdata)%>%
      ggplot(aes(y,Z,col= re[6] ))+geom_point()
    print(p1)

  })

  plot23  <- reactive({
    fdff <- adrarr()$crown_polo
    if (is.null(fdff))
      return(NULL)

    idnum <- fdff[fdff$treeID ==  as.numeric(input$select2),]
    p1<-  ggplot() + geom_sf(data = fdff)+
      geom_sf(data = idnum,col='black')+
      geom_sf_text(data=idnum, aes(label = treeID),col='red',size=6)+
      labs(title = paste0('Total trees count:', length(fdff$treeID))  )+
      theme(plot.title = element_text(size = 20, face = "bold"))
    p1

  })

  plot2  <- reactive({
    fdff <- adrarr()$crown_polo
    if (is.null(fdff))
      return(NULL)

    p1<-  ggplot() + geom_sf(data = fdff)+

      geom_sf_text(data=fdff, aes(label = treeID),col='red',size=1.8)
    p1


  })


  output$predictPlot5  <- renderPlot({
    fdff2 <-  plot23()
    print(fdff2)
  })

  output$downloadsfall  <-  downloadHandler(


    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("spetral_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      withProgress(message = 'Downloading',
                   detail = 'please wait...', value = 0, {
                     png(file, width = input$width_png3, height = input$height_png3, res = input$resolution_PNG3)
                     print(plot2())
                     dev.off()
                   })
    },

    contentType = "application/png"


  )



  rgbplotwithid  <- reactive({
    sfff1 <- adrarr()$sf
    dsf21 <-  getData()
    if (is.null(sfff1))
      return(NULL)
    if (is.null(dsf21))
      return(NULL)
    se2 <-  c(dsf21)

    if(terra::nlyr(se2) < 3){

      tryCatch({
        sp:: plot(se2[[1]] ,col= viridis(200) )
        idnum <- sfff1[sfff1$treeID ==  as.numeric(input$select2),]
        sp::plot(sfff1,border='black', add=T,col=NA,alpha=0.4)
        idnum <- sfff1[sfff1$treeID == as.numeric(input$select2),]
        sp::plot(idnum, add=T,alpha=0.4,col='orange')
        text(idnum, paste(idnum$treeID ),
             cex=1,col='blue' )
      },
      error=function(cond) {
        message("Warning: please upload raster images" )
        print( "Warning: please upload raster images" )
      })
    } else{
      tryCatch({
        library(RStoolbox)
        idnum <- sfff1[sfff1$treeID == as.numeric(input$select2),]
        p <-ggRGB(se2,  stretch = "hist")+
          geom_sf(data = sfff1, fill=NA,col='red' )+
          geom_sf(data = idnum, fill='orange')+
          ggrepel::geom_label_repel(
            data = idnum,
            aes(label = treeID, geometry = geometry),
            stat = "sf_coordinates",
            min.segment.length = 0,
            colour = "red",
            segment.colour = "orange"
          )


        print(p)

      },
      error=function(cond) {
        message("Warning: please upload raster images" )
        print( "Warning: please upload raster images" )
      })
    }

  })


  output$warnpredict  <- renderPrint({
    sfff1 <- adrarr()$sf
    select23 <-  getData()
    rasterDF <-  try (select23   )

    if (inherits(rasterDF, "try-error")){stop(print( "Warning: please upload raster images" ))}


  })

  output$downloadrgball  <- downloadHandler(


    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("spetral_", gsub("/", "-", x), ".png", sep = "")
    } ,
    content = function(file) {
      withProgress(message = 'Downloading',
                   detail = 'please wait...', value = 0, {
                     png(file, width = input$width_png3, height = input$height_png3, res = input$resolution_PNG3)
                     print(rgbplotwithid())
                     dev.off()
                   })
    } ,

    contentType = "application/png"

  )

  output$predictPlo  <- renderPlot({
    rgbplotwithidw<- rgbplotwithid()
    print(rgbplotwithidw)

  })



  finaldata  <- reactive({
    dsf1 <-  data_ext2()
    if (is.null(dsf1))
      return(NULL)
    names(dsf1) <-  gsub('[.]|tif','',names(dsf1))
    dsf1 <- dsf1 %>% dplyr::select( x,y, treeID,Z,area, everything())
    dsf1 <- dsf1 %>% dplyr::select(-ID)  %>% drop_na()
  })

  output$tife_data <- downloadHandler(

    filename = "finaloutput data.rds",
    content = function(file) {
      withProgress(message = 'Downloading',
                   detail = 'please wait...', value = 0, {
                     readr::write_rds(finaldata(), file)
                   })
    }

  )

  output$predictPlot3  <- renderPlot({
    library(tidyverse)
    library(raster)
    library(EBImage)
    library(tools)
    nir <- filter(finaldata(), treeID == as.numeric(input$select2) )
    nir3 <- nir[!names(nir) %in% c('x','y','treeID', 'Z', 'area','ID')]
    library(raster)
    library(RStoolbox)
    library(rasterVis)
    library(viridis)
    nir2 <- sapply(nir3, function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))

    matou_vis2 <- cbind.data.frame(nir[,c('x','y', 'Z', 'area')], nir2)
    matou_vis2 <- matou_vis2 %>% dplyr:: filter(Z > input$heightdata)
    matou_vis2 <- matou_vis2[,colSums(is.na(matou_vis2))<nrow(matou_vis2)]
    t <- try(terra::rast(matou_vis2))

    if (!inherits(t, "try-error")){
      sp:: plot(t ,col= viridis(200))

    }  else {

      e <- extent(matou_vis2[,1:2])
      r <- raster(e, ncol=40, nrow=40)
      x <- rasterize(matou_vis2[, 1:2], r, matou_vis2[,-c(1:2)])
      sp:: plot(x,col= viridis(200))

    }

  })


  output$predictPlot4  <- renderPlot({
    library(tidyverse)
    library(raster)
    library(EBImage)
    library(tools)
    nir <- dplyr::filter(finaldata(), treeID == as.numeric(input$select2) )
    nir3 <- nir[,-c(1:5)]
    library(raster)
    nir2 <- sapply(nir3, function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))

    matou_vis2 <- cbind.data.frame(nir[,c(1:2,4,5)], nir2)
    matou_vis2 <- matou_vis2 %>% dplyr::  filter(Z > input$heightdata)

    rasterDF <- try(terra::rast(matou_vis2))

    if (!inherits(rasterDF, "try-error")){

      if (terra::nlyr(rasterDF) <= 3) {
        print("at least 3 layers needed for RGB plot")
        sp:: plot(rasterDF ,col= viridis(200))

      } else{

        if (terra::nlyr(rasterDF)>3|terra::nlyr(rasterDF)< 5) {
          df <-  ggRGB(rasterDF,  3, 5, 4,
                       stretch = 'lin') + ggtitle(paste0('tree ID-', as.numeric(input$select2)))
          print(df)

        }else{

          df <-  ggRGB(rasterDF, 6, 3, 5,
                       stretch = 'lin') + ggtitle(paste0('tree ID-', as.numeric(input$select2)))
          print(df)


        }
      }
    }  else {
      matou_vis2 <- matou_vis2[,colSums(is.na(matou_vis2))<nrow(matou_vis2)]
      e <- extent(matou_vis2[,1:2])
      r <- raster(e, ncol=40, nrow=40)
      rasterDF <- rasterize(matou_vis2[, 1:2], r, matou_vis2[,-c(1:2)])
      if (terra::nlyr(rasterDF) <= 3) {
        print("at least 3 layers needed for RGB plot")
        sp:: plot(rasterDF ,col= viridis(200))

      } else{

        if (terra::nlyr(rasterDF)>3|terra::nlyr(rasterDF)< 5) {
          df <-  ggRGB(rasterDF,  3, 5, 4,
                       stretch = 'lin') + ggtitle(paste0('tree ID-', as.numeric(input$select2)))
          print(df)

        }else{

          df <-  ggRGB(rasterDF, 6, 3, 5,
                       stretch = 'lin') + ggtitle(paste0('tree ID-', as.numeric(input$select2)))
          print(df)


        }
      }
    }

    library(RStoolbox)
    library(rasterVis)
    library(viridis)

  })

  dfpoint   <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return('please upload las cloud data')
    las_12 <- lapply(inFile$datapath,function(m){
      fdd <- lidR::readLAS(m )
    } )

    las_12

  })

  randomVals <- eventReactive(input$dodo, {
    withProgress(message = 'Ploting',
                 detail = 'May take a while...', value = 0, {
                   sele <- dfpoint()
                   sele <-  sele  %>%  do.call(c,.)
                   lapy <- lapply(sele, function(x){
                     sp::plot(x, bg = "white",size = input$poinsize  , axis = TRUE, legend = TRUE)
                   })
                   lapy

                 })
  })

  output$contents22 <- renderPrint({
    library("lidR")
    library("rgdal")
    library(raster)
    library(tidyverse)
    print("plot with RGL device")
    rma <-  randomVals()
    if (is.null(rma))
      return(NULL)
    print(rma)

    # lapy <- lapply(rma, function(x){
    #   sp::plot(x, bg = "white", , axis = TRUE, legend = TRUE)
    # })
    # lapy
  })







  output$su2  <- renderPrint({

    dsf331 <- finaldata()
    print(dsf331)

  })



  red2 <-  reactive({
    library(raster)
    sele1 <- input$red1
    dsf1 <- raster::raster(sele1$datapath)
  })
  green2 <-  reactive({
    library(raster)
    sele2 <- input$gree1
    dsf2 <- raster::raster(sele2$datapath)
  })
  blue2 <-  reactive({
    library(raster)
    sele3 <- input$blue1
    dsf3 <- raster::raster(sele3$datapath)
  })
  redege2 <-  reactive({
    library(raster)
    sele4 <- input$redege1
    dsf4 <- raster::raster(sele4$datapath)
  })
  NIR2 <-  reactive({
    library(raster)
    sele5 <- input$NIR1
    dsf5 <- raster::raster(sele5$datapath)
  })

  ndviind <- reactive({
    st1 <-  ((NIR2() - red2()) / (NIR2() + red2()))
  })
  osavi <- reactive({
    osavi = ((NIR2() - red2()) * (1 + 0.16)) / (NIR2() + red2() + 0.16)
  })
  gndvi <- reactive({
    gndvi = (NIR2() - green2()) / (NIR2() + green2())
  })
  savi <- reactive({
    savi = ((NIR2() - red2()) * (1 + 0.5)) / ((NIR2() + red2() + 0.5))
  })
  msavi <- reactive({
    msavi = (2 * NIR2() + 1 - sqrt((2 * NIR2() + 1) ^ 2 - 8 * (NIR2() - red2()))) /
      2
  })
  gci <- reactive({
    gci = NIR2() / green2() - 1
  })
  RECI <- reactive({
    RECI = NIR2() / redege2() - 1
  })
  LCI <- reactive({
    LCI = (NIR2() - redege2()) / (NIR2() + red2())
  })
  GRVI <- reactive({
    GRVI = (green2() - red2()) / (green2() + red2())
  })
  MGRVI <- reactive({
    MGRVI = (green2() ^ 2 - red2() ^ 2) / (green2() ^ 2 + red2() ^ 2)
  })
  RGBVI <- reactive({
    RGBVI = (green2() ^ 2 - red2() * blue2()) / (green2() ^ 2 + red2() * blue2())
  })

  NDRE <- reactive({
    NDRE = (NIR2() - redege2()) / (NIR2() + redege2())
  })
  MACI <- reactive({
    MACI = NIR2() / green2()
  })

  ARI <- reactive({
    ARI = green2() / NIR2()
  })

  MARI <- reactive({
    MARI = (green2() ^ (-1) - redege2() ^ (-1)) / NIR2()
  })

  dat243 <- reactive({
    select <- switch(
      input$testh2o,
      "Red" = red2(),
      "Green" = green2(),
      "Blue" = blue2(),
      "Rededage" = redege2(),
      "NIR" = NIR2() ,
      "ndvi" = ndviind(),
      "osavi" = osavi(),
      "gndvi" = gndvi(),
      "savi" = savi(),
      "msavi" = msavi(),
      "gci" = gci(),
      "RECI" = RECI(),
      "LCI" = LCI(),
      "GRVI" = GRVI(),
      "MGRVI" = MGRVI(),
      "NDRE" = NDRE(),
      "MACI" = MACI(),
      "ARI" = ARI(),
      "MARI" = MARI()
    )
  })

  pyt <- reactive({

    st1 <- raster::brick( blue2(), green2(),red2(), redege2(), NIR2())

  })

  output$plotgraph1 <- renderPlot({


    withProgress(message = 'ggRGB Ploting',
                 detail = 'May take a while...', value = 0, {

                   rasterDF <-  try (pyt())
                   if (inherits(rasterDF, "try-error")){stop(print( "Warning: please upload all images" ))}else{
                     library(RStoolbox)
                     library(raster)
                     pyt2 <-  RStoolbox::ggRGB(rasterDF,

                                               stretch  = 'hist')
                     print(pyt2)

                   }
                 })
  })


  output$plotgraph2 <- renderPlot({
    withProgress(message = 'VIs Ploting',
                 detail = 'May take a while...', value = 0, {

                   rasterDF <-  try (dat243())
                   if (inherits(rasterDF, "try-error")){stop(print( "Warning: please upload all images" ))}else{
                     library(RStoolbox)
                     library(raster)
                     sp::plot(dat243())}
                 })


  })

  output$downloadPlotPNG11  <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("spetral_", gsub("/", "-", x), ".tif", sep = "")
    },
    content = function(file) {
      withProgress(message = 'Downloading',
                   detail = 'please wait...', value = 0, {
                     r <- dat243()
                     res <- writeRaster(r, filename=file, format="GTiff", overwrite=TRUE)
                     print(res@file@name)
                   })
    }

  )

  output$downloadPlotPNG22  <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("RGB_", gsub("/", "-", x), ".tif", sep = "")
    },

    content = function(file) {
      withProgress(message = 'Downloading',
                   detail = 'please wait...',
                   value = 0,
                   {
                     r <- pyt()
                     res <- writeRaster(r,
                                        filename = file,
                                        format = "GTiff",
                                        overwrite = TRUE)
                   })
    }
  )


}

app <- shinyApp(ui = ui, server = server)
runApp(app, launch.browser = TRUE)
}


