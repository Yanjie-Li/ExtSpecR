#' 检查和安装所需的 Bioconductor 包
#'
#' @export
check_dependencies <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  pkgs <- c("Biobase", "EBImage")
  if (length(setdiff(pkgs, rownames(installed.packages()))) > 0)
    BiocManager::install(pkgs, ask = FALSE, update = FALSE)

  packages <- c("shinythemes", 'shinyWidgets', 'terra', 'shinyjs', 'RCSF','DT',
                "shinydashboard",'stars','colorspace','readr',
                'sfheaders','sf','exactextractr', 'lidR' ,'tidyverse',
                'viridis', 'rgdal','tictoc',  'ggrepel',
                'raster', 'tools','rasterVis','data.table')

  new.packages <- packages[!(packages %in% utils::installed.packages()[,"Package"])]
  if(length(new.packages)) utils::install.packages(new.packages,quiet =T )

}
