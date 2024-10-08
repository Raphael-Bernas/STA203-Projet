rm(list=ls())
graphics.off()
library(ggplot2)
library(reshape2)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(glmnet)
library(MASS)
library(pls)
library(caret)
library(pROC)
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
theta_hat_test_list <- apply(as.matrix(grid), function(lambda){ridge_regression(lambda)}, MARGIN=1)
coef_test_autre <- coef(glmnet(xtrain_scaled, ytrain_scaled, alpha = 0, lambda = grid, family="gaussian"))
diff_theta_evolution <- apply((as.matrix(coef_test_autre[-1,])-as.matrix(theta_hat_test_list))^2,2,mean)
par(mfrow = c(1, 1))
plot(log(grid), diff_theta_evolution, type = "l", lwd = 4, col = "blue", xlab = "log(lambda)", ylab = "Valeur estimée de la difference L2 des thetas", main = "Variation pour famille gaussienne de la différence L2 des thetas selon lambda ")


set.seed(123)
cv_segments <- cvsegments(n, k = 4)

# Fonction pour calculer le RMSE
rmse_lambda <- function(predictions, actual) {
  return(sqrt(mean((predictions - actual)^2))) 
}

# Initialisation de la matrice pour stocker les RMSE
rmse_cv_lambda_i <- matrix(rep(1, 4 * length(lambda_grid)), nrow = 4)

# Boucle sur les segments de validation croisée
for (i in seq_along(cv_segments)) {
  Y_cv <- ytrain_scaled[cv_segments[[i]]]
  Y_train <- ytrain_scaled[-cv_segments[[i]]]
  X_cv <- xtrain_scaled[cv_segments[[i]], ]
  X_train <- xtrain_scaled[-cv_segments[[i]], ]
  
  # Estimation des coefficients avec régression ridge
  ridge_coef <- coef(glmnet(y = Y_train, x = X_train, family = "gaussian", alpha = 0, lambda = (4/n)*lambda_grid)) # on multiplie lambda par (4/n) qui correspond a la longueur des sous-ensembles, cela est dû a la façon dont lambda est calculé par glmnet 
  teta_cv <- ridge_coef[-1,]  # En variables centrées et réduites, l'intercept est nul
  
  # Prédiction sur les données de validation
  Y_pred_cv <- X_cv %*% teta_cv
  
  # Calcul du RMSE pour chaque valeur de lambda
  rmse_cv_lambda_i[i,] <- apply(Y_pred_cv, MARGIN = 2, FUN = function(predictions) rmse_lambda(predictions, Y_cv))
}

# Calcul de la moyenne des RMSE pour chaque valeur de lambda
rmse_cv_lambda <- apply(rmse_cv_lambda_i, MARGIN = 2, FUN = mean)

# Trouver la valeur optimale de lambda
lambda_opt <- lambda_grid[which.min(rmse_cv_lambda)]

ggplot() + geom_line(mapping = aes(x = log(lambda_grid), y = rmse_cv_lambda)) + 
  xlab("log(lambda)") + ylab("rmse_CV")
lambda_opt

# df contient les degrés de liberté
df <- n - 1
alpha <- 0.05  # Niveau de confiance
CI_down <- sqrt(qchisq(1 - alpha / 2, df) / df) * sqrt(rmse_cv_lambda/n)
CI_up <- sqrt(qchisq(alpha / 2, df) / df) * sqrt(rmse_cv_lambda/n)
data_rmse <- data.frame(lambda = lambda_grid, rmse = rmse_cv_lambda, CI_upper = CI_down, CI_lower = CI_up)

# Tracer le graphique avec les barres d'erreur
ggplot(data_rmse, aes(x = log(lambda), y = rmse)) +
  geom_point() +  # Points pour l'erreur moyenne
  geom_errorbar(aes(ymin = rmse - CI_down, ymax = rmse + CI_up), width = 0.05) +  # Barres d'erreur pour l'intervalle de confiance
  xlab("log(lambda)") + ylab("RMSE") +
  ggtitle("Erreur moyenne avec intervalle de confiance")

fld=rep(0,36)
for (i in 1:36){
  
  fld[i]=(is.element(i,cv_segments$V1))+2*(is.element(i,cv_segments$V2))+3*(is.element(i,cv_segments$V3))+
    4*(is.element(i,cv_segments$V4))
}
cv_model <- cv.glmnet(xtrain_scaled, ytrain_scaled, alpha = 0, lambda = (4/n)*lambda_grid, foldid=fld)

# Extraire les valeurs de lambda et les RMSE moyens à partir du modèle cv.glmnet
lambda_cvglmnet <- cv_model$lambda
rmse_cvglmnet <- cv_model$cvm
ICup_cvglmnet <- cv_model$cvup
IClo_cvglmnet <- cv_model$cvlo

# Créer un dataframe pour les valeurs de lambda, les RMSE et les bornes des intervalles de confiance
data_cvglmnet <- data.frame(lambda = lambda_cvglmnet, rmse = rmse_cvglmnet, CI_upper = ICup_cvglmnet, CI_lower = IClo_cvglmnet)

ggplot(data_cvglmnet, aes(x = log(1/lambda), y = rmse)) +
  geom_point() +  # Points pour l'erreur moyenne
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.05) +  # Barres d'erreur pour l'intervalle de confiance complet
  xlab("log(lambda)") + ylab("RMSE") +
  ggtitle("Erreur moyenne avec intervalle de confiance complet (cv.glmnet)")

lambda_opt # lambda optimale pour validation croisé a la main
ridge_model_final <- glmnet(y = ytrain_scaled, x = xtrain_scaled, alpha = 0, lambda = lambda_opt)
n_test <- length(ytest)
xtest_scaled <- scale(xtest)*(n_test/(n_test-1))
ytest_scaled <- scale(ytest)*(n_test/(n_test-1))
y_pred_test <- predict(ridge_model_final, newx = xtest_scaled)
RMSE_test <- sqrt(mean((ytest_scaled - y_pred_test)^2))
print(paste("Erreur de généralisation (RMSE) CV à la main :", RMSE_test))

lambda_opt_cvglmnet <- lambda_grid[which.min(rmse_cvglmnet)] # lambda optimale trouvé par cvglmnet
lambda_opt_cvglmnet
ridge_model_final <- glmnet(y = ytrain_scaled, x = xtrain_scaled, alpha = 0, lambda = lambda_opt_cvglmnet)
y_pred_test <- predict(ridge_model_final, newx = xtest_scaled)
RMSE_test <- sqrt(mean((ytest_scaled - y_pred_test)^2))
print(paste("Erreur de généralisation (RMSE) cvglmnet :", RMSE_test))

# BONUS : Partie 3
control <- trainControl(method = "cv", number = 5)  # 5-fold cross-validation
# grille des paramètres pour les modèles de régression ridge et lasso
ridge_grid <- expand.grid(alpha = 0, lambda = lambda_grid)
lasso_grid <- expand.grid(alpha = 1, lambda = lambda_grid)

# Ajuster les modèles de régression ridge et lasso
ridge_model <- train(ytrain_scaled ~ ., data = data.frame(xtrain_scaled, ytrain_scaled),
                     method = "glmnet", trControl = control, tuneGrid = ridge_grid)
lasso_model <- train(ytrain_scaled ~ ., data = data.frame(xtrain_scaled, ytrain_scaled),
                     method = "glmnet", trControl = control, tuneGrid = lasso_grid)

# Comparer les performances des modèles
ridge_performance <- ridge_model$results$RMSE
lasso_performance <- lasso_model$results$RMSE

# Afficher les performances
print(paste("Performance du modèle de régression ridge (RMSE) :", mean(ridge_performance)))
print(paste("Performance du modèle de régression lasso (RMSE) :", mean(lasso_performance)))

# La principale différence entre la régression ridge et la régression lasso réside dans la façon dont elles pénalisent les coefficients.
#La régression ridge utilise une pénalité L2 (norme Euclidienne), ce qui conduit à des coefficients réduits vers zéro mais pas exactement à zéro.
#La régression lasso utilise une pénalité L1 (norme de Manhattan), ce qui favorise la parcimonie en forçant de nombreux coefficients à zéro, ce qui permet la sélection de variables.


## Partie 4 : Régression logistique pénalisé

z <- (ytrain > 88) - 0
ztest <- (ytest > 88) - 0
sum(z==1)/sum(z==0)
sum(ztest==1)/sum(ztest==0)
# les jeux sont plutôt équilibré

ridge_cv <- cv.glmnet(x = xtrain_scaled, y = z, family = "binomial", alpha = 0,
                      foldid = fld, type.measure="class")
lasso_cv <- cv.glmnet(x = xtrain_scaled, y = z, family = "binomial", alpha = 1,
                      foldid = fld, type.measure="class")


par(mfrow = c(1, 2))  
plot(ridge_cv, main = "modèle de ridge")
plot(lasso_cv, main = "modèle de lasso")

ridge_prob <- predict(ridge_cv, newx = xtest_scaled, s = "lambda.min", type = "response")
lasso_prob <- predict(lasso_cv, newx = xtest_scaled, s = "lambda.min", type = "response")
# Bayes
ridge_pred <- ifelse(ridge_prob > 0.5, 1, 0)
lasso_pred <- ifelse(lasso_prob > 0.5, 1, 0)
# Erreur de généralisation (taux d'erreur) pour chaque modèle
ridge_error <- mean(ridge_pred != ztest)
lasso_error <- mean(lasso_pred != ztest)
print(paste("Erreur de généralisation pour la régression logistique pénalisée en ridge :", ridge_error))
print(paste("Erreur de généralisation pour la régression logistique pénalisée en lasso :", lasso_error))

# Tracer les courbes roc
par(mfrow=c(1,1))
plot(roc.glmnet(ridge_cv,newx=xtest_scaled,newy=ztest,family="binomial"),type='l', col = "blue")
lines(roc.glmnet(lasso_cv,newx=xtest_scaled,newy=ztest,family="binomial"),type='l', col = "red")
abline(0, 1, col = "black")
legend("bottomright", legend = c("Ridge", "Lasso"), col = c("blue", "red"), lty = 1, cex = 0.8)
# Lasso gagne devant ridge : en effet aire entre droite noir et courbe de lasso plus grande que c'elle pour ridge

