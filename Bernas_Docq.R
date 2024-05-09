rm(list=ls())
graphics.off()
library(ggplot2)
library(reshape2)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(glmnet)
library(MASS)
# Vous pouvez changer la variable path pour changer le working directory
path = "C:/Users/R/OneDrive/Bureau/Travail/Ensta_Paris/Math/STA203/Projet-STA203"
setwd(path)

## Partie 2 : Analyse exploratoire

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
acp <- PCA(xtrain, ncp = 5) #  graph = FALSE

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
  axes <- res$ind$coord[,1:nr]
  vars <- t(res$var$coord[,1:nr])
  if(nr!=1){
    eigenvalue_inv<- diag(1/sqrt(res$eig[1:nr,1]))
    reconstructed <- axes%*%eigenvalue_inv%*%vars
  }
  else{
    eigenvalue_inv<- 1/sqrt(res$eig[1:nr,1])
    reconstructed <- eigenvalue_inv*axes%*%vars
  }
  
  # Réduire les axes selon les écarts-types et les moyennes des variables explicatives
  cloud <- t(t(reconstructed) * Xsd + Xm)
  
  # Retourner les axes reconstruits
  return(cloud)
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
reconstruction_totale <- reconstruct(res = acp, nr = 1, Xm = colMeans(xtrain), Xsd = apply(xtrain, 2, sd))
par(mfrow = c(1, 1)) 
matplot(t(reconstruction_totale), type = 'l', main = "Courbes des Spectres",
        xlab = "Fréquences", ylab = "Mesures")
rmse_total <- rmse(xtrain, reconstruction_totale)
mae_total <- mae(xtrain, reconstruction_totale)
cat("RMSE de la reconstruction totale :", rmse_total, "\n")
cat("MAE de la reconstruction totale :", mae_total, "\n\n")

# Tracer les reconstructions pour chaque nr
par(mfrow = c(2, 3))  # Diviser l'espace graphique en 2 lignes et 3 colonnes
for (nr in c(1, 2, 3, 4, 5)) {
  # Calcul de la reconstruction pour nr
  reconstruction <- reconstruct(res = acp, nr = nr, Xm = colMeans(xtrain), Xsd = apply(xtrain, 2, sd))
  
  # Calcul de l'erreur
  rmse_val <- rmse(xtrain, reconstruction)
  mae_val <- mae(xtrain, reconstruction)
  
  # Tracer la reconstruction avec le titre
  plot(reconstruction, main = paste("Reconstruction pour nr =", nr, "\nRMSE =", round(rmse_val, 2), "MAE =", round(mae_val, 2)))
}
par(mfrow = c(2, 3))  # Diviser l'espace graphique en 2 lignes et 3 colonnes
for (nr in c(1, 2, 3, 4, 5)) {
  # Calcul de la reconstruction pour nr
  reconstruction <- reconstruct(res = acp, nr = nr, Xm = colMeans(xtrain), Xsd = apply(xtrain, 2, sd))
  
  # Calcul de l'erreur
  rmse_val <- rmse(xtrain, reconstruction)
  mae_val <- mae(xtrain, reconstruction)
  
  # Tracer la reconstruction avec le titre
  matplot(t(reconstruction), type = 'l', main = paste("Spectre pour nr =", nr, "\nRMSE =", round(rmse_val, 2), "MAE =", round(mae_val, 2)), xlab = "Fréquences", ylab = "Mesures")
}

## Partie 3 : Régression pénalisé

grid=10^seq(6,-10,length=100)
# Estimer le modèle de régression ridge
ridge_model <- glmnet(xtrain, ytrain, alpha = 0, lambda = grid)

# Extraire les coefficients du modèle
coefficients <- coef(ridge_model)

# Extraire les valeurs estimées du paramètre d'intercept
intercept_values <- coefficients[1, ]

# Tracer la variation de la valeur estimée du paramètre d'intercept
par(mfrow = c(1, 1))
plot(log(grid), intercept_values, type = "l", lwd = 4, col = "blue", xlab = "log(lambda)", ylab = "Valeur estimée de l'intercept", main = "Variation de l'intercept avec lambda (ridge)")

# Intercept calculer avec les autres estimations
Intercept_estime=apply(-(as.matrix(xtrain)%*%coef(ridge_model)[-1,]-as.vector(ytrain)),mean,MARGIN=2)
lines(log(grid), Intercept_estime, type = "l", lwd =2, col = "red", xlab = "log(lambda)", ylab = "Valeur estimée de l'intercept", main = "Variation de l'intercept avec lambda (ridge)")
legend("topright", legend = c("Intercepte du modèle", "Intercepte recalculé"), lty = 1, col = c("blue", "red"), lwd = c(4, 2))

# On centre :
xtrain_tilde = scale(xtrain, scale=FALSE)
ridge_model_tilde1 <- glmnet(xtrain_tilde, ytrain, alpha = 0, lambda = grid)
coefficients_tilde1 <- coef(ridge_model_tilde1)
intercept_values_tilde1 <- coefficients_tilde1[1, ]
par(mfrow = c(1, 1))
plot(log(grid), intercept_values_tilde1, type = "l", lwd = 3, col = "blue", xlab = "log(lambda)", ylab = "Valeur estimée de l'intercept", main = "Variation de l'intercept avec lambda (ridge)", ylim = c(-10, 100))

ytrain_tilde = scale(ytrain, scale=FALSE)
ridge_model_tilde2 <- glmnet(xtrain_tilde, ytrain_tilde, alpha = 0, lambda = grid)
coefficients_tilde2 <- coef(ridge_model_tilde2)
intercept_values_tilde2 <- coefficients_tilde2[1, ]
lines(log(grid), intercept_values_tilde2, type = "l", lwd = 3, col = "red", xlab = "log(lambda)", ylab = "Valeur estimée de l'intercept", main = "Variation de l'intercept avec lambda (ridge)")
legend("topright", legend = c("Intercepte du modèle X centré", "Intercepte du modèle X et Y centré"), lty = 1, col = c("blue", "red"), lwd = c(4, 2))

#On centre-réduit
n = length(xtrain[,1])
xtrain_scaled = scale(xtrain)*sqrt(n/(n-1))
ytrain_scaled = scale(ytrain)*sqrt(n/(n-1))
svd_xtrain <- svd(xtrain_scaled)
u <- svd_xtrain$u
v <- svd_xtrain$v
s <- svd_xtrain$d
A0 <- v %*% diag((1 / s)) %*% t(u)
mae(A0 %*% xtrain_scaled, diag(rep(1,401))) # MAE de I - A0*X : proche de 0 donc on à presqu'une inverse

theta = A0%*%ytrain_scaled # valeur observé de theta
head(theta)
head(A0)

ridge_model_scaled <- glmnet(xtrain_scaled, ytrain_scaled, alpha = 0, lambda = grid)
coefficients_scaled <- coef(ridge_model_scaled)

bias_ridge_estimate <- apply((as.matrix(coefficients_scaled[-1,])-as.vector(theta)),2,mean)
par(mfrow = c(1, 1))
plot(log(grid), bias_ridge_estimate, type = "l", lwd = 4, col = "blue", xlab = "log(lambda)", ylab = "Valeur estimée du biais de ridge", main = "Variation du biais éstimé avec lambda ")

# lm.ridge
res.ridge <- lm.ridge(ytrain_scaled ~ xtrain_scaled, lambda = grid)
coefficients_lmridge <- res.ridge$coef
# On calcul ici un RSE entre les prédictions du theta obtenu a la question précédente qui est une estimation correcte de theta quand lambda tend vers 0 et les prédictions des thetas de lm.ridge
MSE_ridge_estimate <- apply((as.matrix((xtrain_scaled %*% coefficients_lmridge)) - as.vector(xtrain_scaled %*% theta))^2, 2, mean)
par(mfrow = c(1, 1))
plot(log(grid), MSE_ridge_estimate, type = "l", lwd = 4, col = "blue", xlab = "log(lambda)", ylab = "Valeur estimée du RSE de ridge", main = "Variation du RSE éstimé avec lambda ")

ridge_regression <- function(lambda) {
  XTX <- t(xtrain_scaled) %*% xtrain_scaled
  XTX_lambda <- XTX + lambda * diag(ncol(xtrain_scaled))
  XTY <- t(xtrain_scaled) %*% (ytrain_scaled - mean(ytrain_scaled)*rep(1,length(ytrain_scaled)))
  theta_hat <- solve(XTX_lambda) %*% XTY
  return(theta_hat)
}
lambda <- 0.001
theta_hat_list <- apply(as.matrix(grid), function(lambda){ridge_regression(lambda/2)}, MARGIN=1)
head(theta_hat_list[,1])
head(coefficients_scaled[-1,1])
diff_theta_evolution <- apply((as.matrix(coefficients_scaled[-1,])-as.matrix(theta_hat_list))^2,2,mean)
par(mfrow = c(1, 1))
plot(log(grid), diff_theta_evolution, type = "l", lwd = 4, col = "blue", xlab = "log(lambda)", ylab = "Valeur estimée de la difference L2 des thetas", main = "Variation de la différence L2 des thetas selon lambda ")
# On a une faible difference donc une méthode calculé très proche de la méthode utilisé par glmnet
