source("functions.R")
source("transfer.R")
library(ggplot2)
library(mdatools)
library(dplyr)
library(rafalib)
library(plotly)
library(pls)
library(robustbase)
library(parallel)
mypar(1,1)
pls.options(parallel=6)

raw_calset <- read.csv("CALSET.csv", header = TRUE)
y_response <- as.numeric(raw_calset[,3])
spectra <- raw_calset

#Remove columns without data
rcol <- c(1,2,3)
spectra <- spectra[,-rcol]
wavelength <- seq(300,1100, length=ncol(spectra))
spectra <- apply(spectra, 2, as.numeric)

#############Data preprocessing###############################
#Ordering data
ordered_i <- order(y_response)
spectra <- spectra[ordered_i,]
y_response <- y_response[ordered_i]


#Select wavelengths
#intervals <- ipls(spectra, y_response, glob.ncomp=7, int.num=32, method="forward")
#summary(intervals)
#plot(intervals)
#spectra <- chop_spectra(spectra,wavelength,750,1000)
#Load from matlab FIXME
#spectra <- read.csv("750-1000_outliers.csv", header=FALSE)
#spectra <- apply(spectra, 2, as.numeric)
#y_response <- read.csv("y_response_outliers.csv", header=TRUE)
#y_response <- y_response[,1]

wl1 <- 750
wl2 <- 950

#Pre-processing
spectra <- snv(spectra)
spectra <- my_savgol(spectra, 21, 2, 2)
spectra <- chop_spectra(spectra,wavelength,wl1,wl2)

#spectra_plot_y(spectra, wavelength, y_response)
#Select wavelengths
#summary(intervals)
#Visualize additive and multiplicative effects
#means <- apply(spectra, 2, mean)
#spectra_plot_y(spectra, wavelength, y_response)

#############End of data preprocessing########################

#############Spectra visualization############################
#spectra_plot_y(spectra, wavelength, y_response)
#############END Spectra visualization############################

#####------------PCA------.#######
#spectra_pca <- prcomp(spectra, center=FALSE, scale.= FALSE)

#pca_data <- data.frame(PC1=spectra_pca$x[,1], PC2=spectra_pca$x[,2],
#                       PC3=spectra_pca$x[,3],PC4=spectra_pca$x[,4], Y=y_response)
#scores plot
#scores <- ggplot() + 
#        geom_point(data=pca_data, aes(x=PC1, y=PC2, colour=Y)) + 
#        scale_color_gradientn(colors=rainbow(5))
#plot(scores)

#plot_ly(x=pca_data$PC1, y=pca_data$PC2, z=pca_data$PC3, type="scatter3d", mode="markers", color=pca_data$Y)

#summary(spectra_pca)
#Screeplot
#screeplot(spectra_pca, main="Screeplot", npcs=10, type="lines")
######END of PCA################################################

ncomp <- 10

spectra_data <- data.frame(y_response=y_response, NIR=I(spectra))

pls_model <- plsr(y_response ~ NIR, data=spectra_data, ncomp=ncomp) 
pls_cv <- crossval(pls_model, segments=5)
plot(RMSEP(pls_cv), legendpos="topright")
#summary(pls_cv, what="validation")

scores_data <- data.frame(LV1=pls_model$scores[,2],LV2=pls_model$scores[,3], LV3=pls_model$scores[,3],Y=y_response)
scores_data <- data.frame(LV=I(pls_model$scores), Y=y_response)
scores <- ggplot() + 
        geom_point(data=scores_data, aes(x=LV[,1], y=LV[,2], colour=Y)) + 
        scale_color_gradientn(colors=rainbow(5))
#plot(scores)

#plot_ly(x=pls_model$scores[,1], y=pls_model$scores[,2], z=pls_model$scores[,3], type="scatter3d", mode="markers", color=y_response)

#plot(pls_model$residuals[,1,ncomp])
RMSEP <- sqrt(pls_model$residuals[,1,ncomp]%*%pls_model$residuals[,1,ncomp]/nrow(spectra))
print(RMSEP)

validation <- generate_validation(msc_reference)
#spectra_plot(validation, wavelength)
validation <- chop_spectra(validation,wavelength,wl1,wl2)

predictions <- predict(pls_model, ncomp=ncomp, newdata=validation)

results <- read.csv("Ytest.csv", header=TRUE)
results <- results[,1]

plot(results,predictions)

RMSEP <- sqrt((predictions-results)%*%(predictions-results)/length(predictions))
print(RMSEP)
write.csv(predictions, "predictions.csv")
