#***************************************************
### Packages  ----
#***************************************************

if (!require("tcltk2")) {
  install.packages("tcltk2", dependencies = TRUE)
#  library(tcltk2)
}
if (!require("flowCore")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("flowCore")
#  library(flowCore)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
#  library(ggplot2)
}
if (!require("cowplot")) {
  install.packages("cowplot", dependencies = TRUE)
#  library(cowplot)
}
if (!require("scales")) {
  install.packages("scales", dependencies = TRUE)
#  library(scales)
}
if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
#  library(tidyverse)
}
if (!require("ggjoy")) {
  install.packages("ggjoy", dependencies = TRUE)
#  library(ggjoy)
}
if (!require("binovisualfields")) {
  install.packages("binovisualfields", dependencies = TRUE)
#  library(binovisualfields)
}
if(!require(ggpubr)){ 
  install.packages("ggpubr", dependencies = TRUE)
#  library(ggpubr)
}
if(!require(pastecs)){ 
  install.packages("pastecs", dependencies = TRUE)
#  library(pastecs)
}


