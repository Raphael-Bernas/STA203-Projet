rm(list=ls())
graphics.off()
library(ggplot2)
library(reshape2)
library(corrplot)
library(FactoMineR)
library(factoextra)
# Vous pouvez changer la variable path pour changer le working directory
path = "C:/Users/R/OneDrive/Bureau/Travail/Ensta_Paris/Math/STA203/Projet-STA203"
setwd(path)

# Récupération des jeux de données
df_test=read.table("gasolineTest.txt",header=T)
dim(df_test)
df_train=read.table("gasolineTrain.txt",header=T)
dim(df_train)

xtrain <- as.matrix(df_train[, -1])  # Exclure la première colonne (indice d'octane)
xtest <- as.matrix(df_test[, -1])    # Exclure la première colonne (indice d'octane)

ytrain <- df_train$octane
ytest <- df_test$octane

par(mfrow = c(4, 6))  

for (i in 1:24) {
  boxplot(xtrain[, i], main = colnames(xtrain)[i], ylab = "Valeurs")
}

# Courbes des spectres
par(mfrow = c(1, 1)) 
matplot(t(xtrain), type = 'l', main = "Courbes des Spectres",
        xlab = "Fréquences", ylab = "Mesures")

# Corrélation 
correlation_matrix <- cor(xtrain)
heatmap(correlation_matrix, 
        main = "Corrélation entre les Mesures aux Différentes Fréquences",
        xlab = "Fréquences", ylab = "Fréquences")
corrplot(correlation_matrix,method = "circle") # Attention c'est long a calculé

# ACP
acp <- PCA(xtrain) #  graph = FALSE

# Graphe des valeurs propres
par(mfrow = c(1, 1))
eigenvalues <- acp$eig[, "eigenvalue"]
barplot(acp$eig[, "percentage of variance"], main = "Contribution des Valeurs Propres",
        xlab = "Composantes Principales", ylab = "Pourcentage de Variance Expliqué")

barplot(acp$eig[, "cumulative percentage of variance"], main = "Pourcentage Cumulatif de Variance",
        xlab = "Composantes Principales", ylab = "Pourcentage Cumulatif de Variance Expliqué")

# Plot des nuages dans les six premiers axes principaux
# Définir les combinaisons d'axes à afficher
combinaisons_axes <- list(c(1, 2), c(3, 4), c(5, 1))

# Tracer les nuages de points pour chaque combinaison d'axes
for (comb in combinaisons_axes) {
  plot <- fviz_pca_ind(acp, axes = comb, geom = "point")
  print(plot)
}

reconstruct <- function(res, nr, Xm, Xsd) {
  # Sélectionner les premiers nr axes de l'ACP
  axes <- res$ind$coord[, 1:nr]
  
  # Réduire les axes selon les écarts-types et les moyennes des variables explicatives
  axes_reduced <- t(apply(axes, 1, function(x) (x * Xsd) + Xm))
  
  # Retourner les axes reconstruits
  return(axes_reduced)
}
# Fonction pour calculer l'erreur quadratique moyenne (RMSE)
rmse <- function(x, y) {
  return(sqrt(mean((x - y)^2)))
}

# Fonction pour calculer l'erreur en valeur absolue (MAE)
mae <- function(x, y) {
  return(mean(abs(x - y)))
}

# Vérification de la reconstruction totale du nuage avec xtrain
reconstruction_totale <- reconstruct(res = acp, nr = ncol(xtrain), Xm = colMeans(xtrain), Xsd = apply(xtrain, 2, sd))
rmse_total <- rmse(xtrain, reconstruction_totale)
mae_total <- mae(xtrain, reconstruction_totale)
cat("RMSE de la reconstruction totale :", rmse_total, "\n")
cat("MAE de la reconstruction totale :", mae_total, "\n\n")

# Tracer les reconstructions pour chaque nr
par(mfrow = c(2, 3))  # Diviser l'espace graphique en 2 lignes et 3 colonnes
for (nr in c(1, 2, 3, 4, 5, 39)) {
  # Calcul de la reconstruction pour nr
  reconstruction <- reconstruct(res = acp, nr = nr, Xm = colMeans(xtrain), Xsd = apply(xtrain, 2, sd))
  
  # Calcul de l'erreur
  rmse_val <- rmse(xtrain, reconstruction)
  mae_val <- mae(xtrain, reconstruction)
  
  # Tracer la reconstruction avec le titre
  plot(reconstruction, main = paste("Reconstruction pour nr =", nr, "\nRMSE =", round(rmse_val, 2), "MAE =", round(mae_val, 2)))
}
