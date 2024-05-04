rm(list=ls())
graphics.off()
# Vous pouvez changer la variable path pour changer le working directory
path = "C:/Users/R/OneDrive/Bureau/Travail/Ensta_Paris/Math/STA203/Projet-STA203"
setwd(path)

# Récupération des jeux de données
df_test=read.table("gasolineTest.txt",header=T)
dim(df_test)
df_train=read.table("gasolineTrain.txt",header=T)
dim(df_train)
