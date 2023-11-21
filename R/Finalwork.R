library(tidyverse)
library(dplyr)
library(caret)
library(class)
library(gmodels)
library(psych)
devtools::install_github('AndresMCB/DynamicCancerDriverKM')
library(DynamicCancerDriverKM)

datanormal <-(DynamicCancerDriverKM::BRCA_normal)
datapt <-(DynamicCancerDriverKM::BRCA_PT)

final_data <- bind_rows(datanormal, datapt)

