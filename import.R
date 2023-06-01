## Importar paquetes

packages <- c("plyr", "xlsx", "FactoMineR", "factoextra", "qqman", "devtools", "tidyverse", "ggVennDiagram", "remotes")

install.packages(setdiff(packages, rownames(installed.packages())))  

if(system.file(package='fastman')==''){
    devtools::install_github('kaustubhad/fastman',build_vignettes = TRUE)
}

if(system.file(package='summarytools')==''){
    remotes::install_github("dcomtois/summarytools", ref="dev-current")
}
