library(tidyverse)
library(dplyr)
library(caret)
library(class)
library(gmodels)
library(psych)
library(tidyr)
devtools::install_github('AndresMCB/DynamicCancerDriverKM')
devtools::install_github('AndresMCB/AMCBGeneUtils')
library(DynamicCancerDriverKM)
#step 1
datanormal <-(DynamicCancerDriverKM::BRCA_normal)
datapt <-(DynamicCancerDriverKM::BRCA_PT)
#Step 2
final_data <- bind_rows(datanormal, datapt)
#Step 3
porcentaje_menor_10 <- colMeans(final_data < 700, na.rm = TRUE)
columnas_a_eliminar <- names(porcentaje_menor_10[porcentaje_menor_10 >= 0.8])
final_data_filtrado <- final_data[, !names(final_data) %in% columnas_a_eliminar]
final_data_filtrado2 <- final_data
#Step 4
PPI<-(DynamicCancerDriverKM::PPI)
Data_PPI<- PPI %>%
  pivot_longer(cols = c(`Input-node Gene Symbol`, `Output-node Gene Symbol`), names_to = "variable", values_to = "gen") %>%
  group_by(gen, variable) %>%
  summarise(frecuencia = n()) %>%
  pivot_wider(names_from = variable, values_from = frecuencia, values_fill = 0)
Data_PI<- Data_PPI %>%
  mutate(total_mode = `Input-node Gene Symbol`, `Output-node Gene Symbol`) %>%
  select(total_mode) %>%
  arrange(desc(total_mode))
print(Data_PI)
PPI_INIC_FINAL<-colnames(final_data_filtrado)[ 8:ncol(final_data_filtrado)]
aux2 <- AMCBGeneUtils::changeGeneId(PPI_INIC_FINAL, from = "Ensembl.ID")
names(final_data_filtrado)[8:11357] <- aux2$HGNC.symbol
FINAL_GENES<- colnames(final_data_filtrado)
Data_PI_Filtrado<- Data_PI%>%
filter(gen %in% FINAL_GENES)
#Step 5
#PART A
#KNN MODEL




