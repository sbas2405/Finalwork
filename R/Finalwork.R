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

porcentaje_menor_10 <- colMeans(final_data < 700, na.rm = TRUE)
columnas_a_eliminar <- names(porcentaje_menor_10[porcentaje_menor_10 >= 0.8])
final_data_filtrado <- final_data[, !names(final_data) %in% columnas_a_eliminar]
final_data_filtrado2 <- final_data

PPI<-(DynamicCancerDriverKM::PPI)


