library(car)
library(corrplot)
library(dplyr)
library(rafalib)

set.seed(1)
mypar(2,2)
Rcrit <- function(probability,n)
{
        df <- n-2
        rcrit <- qt(1-probability/2,df)/sqrt(qt(1-probability/2,df)^2+df) 
        
        return(rcrit)
}

#################################################3#
data <- as.data.frame(read.csv("data.csv"))

#Removing text
data <- data %>% select(-Solvent) %>% select(-ID)

#Model from article
modelX <- lm(b.p. ~ HBA + HBD + RB, data=data) 
#Checking for outliers
#boxplot(data$b.p., main="Boiling point of glycerol based solvents", ylab ="Boiling point (ºC)")

# Checking for correlation
corr <- round(cor(data),5)
#corrplot(corr, type = "lower", order = "original", method="number", tl.col = "black", tl.srt = 35, number.cex=0.75)
highcorr <- c()
cutoff <- 0.95

#Deleting highly correlated data
for(i in 1:(ncol(corr)-2))
    for(j in (i+1):(nrow(corr)-1))
            if(abs(corr[i,j]) >= cutoff)
            {
                    if(abs(corr[nrow(corr),j]) >= abs(corr[nrow(corr),i]))
                            highcorr <- append(highcorr,i)
                    else
                            highcorr <- append(highcorr,j)
            }            

highcorr <- unique(highcorr)
data <- data[,-highcorr]
#corr <- round(cor(data),5)
#corrplot(corr, type = "lower", order = "original", method="number", tl.col = "black", tl.srt = 35, number.cex=0.75)

#Rcrit
#corr <- round(cor(data),5)
#rcrit <- Rcrit(0.05,nrow(data))
#i <- which(abs(corr[,ncol(corr)]) >= rcrit)
#data <- data[,i]

#Dividing data into training and test sets
index <- 1:nrow(data) 
index <- sample(index,length(index))
size <- 0.70

train_index <- index[1:round(length(index)*size)]
test_index <- index[(round(length(index)*size)+1):length(index)]

train <- data[train_index,]
test <- data[test_index,]

#Start modelling
pval <- c()
Fval <- c()
Fpval <- c()
unsig_var <- c()
Ra <- c()
############################
model <- lm(b.p. ~ ., data=train)

Ra <- append(Ra,summary(model)$adj.r.squared)
Fval <- append(Fval, summary(model)$fstatistic[1])
a <- summary(model)
Fpval <- append(Fpval, pf(a$fstatistic[1],a$fstatistic[2],a$fstatistic[3],lower.tail=FALSE))

i <- which.max(summary(model)$coefficients[-1,4])

unsig_var <- append(unsig_var, names(i))
pval <- append(pval, summary(model)$coefficients[i+1,4])

repeat
{
        #Remove variable
        trainprev <- train
        modelprev <- model
        testprev <- test

        train <- train[,-i]
        test <-test[,-i]

        model <- lm(b.p. ~ ., data=train)

        Ra <- append(Ra,summary(model)$adj.r.squared)
        Fval <- append(Fval, summary(model)$fstatistic[1])
        a <- summary(model)
        Fpval <- append(Fpval, pf(a$fstatistic[1],a$fstatistic[2],a$fstatistic[3],lower.tail=FALSE))

        i <- which.max(summary(model)$coefficients[-1,4])

        unsig_var <- append(unsig_var, names(i))
        pval <- append(pval, summary(model)$coefficients[i+1,4])
        #Remove variable
        if((pval[length(pval)] <= 0.05) || (Ra[length(Ra)] <= Ra[(length(Ra)-1)]))
                break
}
#Saving the previous model before deleting the last unsignificant term
train <- trainprev
test <- testprev
model <- modelprev

#VIF
#i <- c(4,5,6,7,8)
i <- c(4,5)
train <- train[,-i]
test <- test[,-i]
model <- lm(b.p. ~ ., data=train)
print(summary(model))
#Looking for outliers in model
std_res <- scale(model$residuals)

plot(y=std_res, x=train$b.p., main="Residual plot for training set (standardized residuals vs. boiling point (ºC))",
     xlab = "Boiling point (ºC)", ylab = "Std. residuals", col = "blue", ylim=c(-4,4))
abline(h=c(0,-2,2,-3,3),col = c("green","orange","orange","red","red"), lty=c(2,2,2))
legend("right", inset=0.01,legend = c("Training set", "2 std.dev", "3 std.dev"), col = c("blue","orange","red"), pch = c(1,NA,NA), lty = c(NA,2,2),bg="white")


#Predictions on test set
predic <- predict.lm(object=model, newdata=test)
std_res_test <- scale((test[,ncol(test)]-predic))

plot(y=std_res_test, x=test$b.p., main="Residual plot for test set (standardized residuals vs. boiling point (ºC))",
     xlab = "Boiling point (ºC)", ylab = "Std. residuals", col = "blue", ylim=c(-4,4))
abline(h=c(0,-2,2,-3,3),col = c("green","orange","orange","red","red"), lty=c(2,2,2))
legend("right", inset=0.01,legend = c("Test set", "2 std.dev", "3 std.dev"), col = c("blue","orange","red"), pch = c(1,NA,NA), lty = c(NA,2,2),bg="white")

#Calculate R2 for train and test sets
test_R2 <- cor(test$b.p.,predic)^2
#test_R2 <- 1 - (1-test_R2)*(nrow(test)-1)/(nrow(test) - (ncol(test)-1) -1) 
train_R2 <- summary(model)$adj.r.squared

num <- predic-test$b.p.
den <- train$b.p. - mean(predict.lm(object=model, newdata=train))
Q2 <- 1-sum(num^2)/sum(den^2)

#Scatterplots
plot(x = fitted(model), y= train$b.p., main="Estimated boiling point vs. experimental boiling point" ,xlab = "Estimated boiling point (ºC)", ylab = "Boiling point (ºC)",col = "blue")
points(x = predic, y = test$b.p., col = "red")
abline(0,1)
legend("bottomright", inset=0.01, legend = c(paste("Training set, R² =",round(train_R2,3)), paste("Test set, R² =",round(test_R2,3)), 
                             paste("Predictive accuracy, Q² =", round(Q2,3))), col = c("blue", "red"),pch =c(1,1,NA))



print(vif(model))
results <- cbind(Ra, Fval, Fpval, unsig_var, pval)
write.csv(results, "results.csv")


#corr <- round(cor(train),2)                                         
#corrplot(corr, type = "lower", order = "original", method="number", tl.col = "black", tl.srt = 35, number.cex=1)

