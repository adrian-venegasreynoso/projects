library(dplyr)
library(rafalib)
library(car)
mypar(2,2)
set.seed(1)

bit_flip <- function(bit)
{
        if(bit == 0)
                bit = 1
        else
                bit = 0
}

#Bit flip
mutate <- function(parent,mutation_pool)
{
        mutant <-vector("numeric", length(parent))
        repeat
        {
                for(i in 1:length(parent))
                {
                        if(sample(mutation_pool,1) == 1)
                                mutant[i] <- bit_flip(parent[i])
                        else
                                mutant[i] <- parent[i]
                }
                if(sum(mutant) > 0)
                        break
        }

        return(mutant)
}

generate_chromosome <- function(no_genes)
{
        repeat{      
                chromosome <- sample(bits,no_genes, replace=TRUE)
                if(sum(chromosome > 0))
                        break
        }
        return(chromosome)
}

crossover <- function(parent1, parent2)
{
        crossover_child <- vector("numeric", length(parent1))
        splice_index <- 2:(length(parent1)-1) #Don't include element 1 or last 
        repeat
        {
                crossover_gene <- sample(splice_index,1)
                crossover_child <- append(parent1[1:crossover_gene],parent2[(crossover_gene+1):length(parent2)])
                if(sum(crossover_child) > 0)
                        break
        }
   
        return(crossover_child)
}


#Generationg vector for mutation chance
mutation_chance <-  1 #In percent
mutation_1s <- replicate(mutation_chance,1) 
mutation_0s <- vector("numeric", (100-mutation_chance))
mutation_pool <- append(mutation_0s, mutation_1s)

bits <- c(0,1)

########### Linear regression ##################
data <- as.data.frame(read.csv("data.csv"))                                                                                                                                                                       
#Removing text
data <- data %>% select(-Solvent) %>% select(-ID)

# Checking for correlation
corr <- round(cor(data),5)
highcorr <- c()
cutoff <- 0.99
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

#Dividing data into training and test sets
index <- 1:nrow(data) 
index <- sample(index,length(index))
size <- 0.70

train_index <- index[1:round(length(index)*size)]
test_index <- index[(round(length(index)*size)+1):length(index)]

train <- data[train_index,]
test <- data[test_index,]

#############Implementation
no_genes <- ncol(train)-1
no_mutants <- 6
no_crossover <- 2
no_elitist <- 2
F <- c()
elitist_i <- vector("numeric",no_elitist)

for(i in 1:200)
{
        if(i == 1)
        {
                generation_size <- 25
                chromosome_matrix <- matrix(,no_genes,generation_size)
                f <- vector("numeric", generation_size)

                for(i in 1:generation_size)
                        chromosome_matrix[,i] <- generate_chromosome(no_genes)
        }
        else
        {
                generation_size <- no_mutants + no_crossover + no_elitist
                chromosome_matrix <- matrix(,no_genes,generation_size)
                f <- vector("numeric", generation_size)

                for(k in 1:generation_size)
                        if(k <= no_elitist)
                                chromosome_matrix[,k] <- chromosome_matrix_prev[,elitist_i[k]]
                        else if (k <= (no_elitist + no_crossover))
                                 {
                                         if(k == (no_elitist+1))
                                            {
                                                parent1 <- chromosome_matrix_prev[,elitist_i[1]]
                                                parent2 <- chromosome_matrix_prev[,elitist_i[2]]
                                                chromosome_matrix[,k] <- crossover(parent1, parent2)
                                           }
                                         else
                                         {
                                                parents_i <- sample(generation_size,2)
                                                parent1 <- chromosome_matrix_prev[,parents_i[1]]
                                                parent2 <- chromosome_matrix_prev[,parents_i[2]]
                                                chromosome_matrix[,k] <- crossover(parent1, parent2)
                                         }


                                 }
                        else
                        {
                                mutant_i <- sample(generation_size,1)
                                chromosome_matrix[,k] <- mutate(chromosome_matrix_prev[,mutant_i],mutation_pool)
                        }

        }
        #Test the models
        for(j in 1:generation_size)
        {
                expressed_genes <- which(chromosome_matrix[,j] == 1)
                expressed_genes <- append(expressed_genes, ncol(train)) #Adding the Y

                chromosome_train <- train[,expressed_genes]
                 
                model <- lm(b.p. ~ ., data=chromosome_train)

                if(length(expressed_genes) == 2)
                        f[j] <- -99999999999999999999 
                else    
                        f[j] <- summary(model)$adj.r.squared -mean(vif(model))#- (length(expressed_genes)-1)/(ncol(data)-1)/100 - (mean(vif(model))-1)
        }

        F <- append(F, max(f))

        sorted <- sort(f, decreasing=TRUE)

        for(j in 1:no_elitist)
                elitist_i[j] <- which(f == sorted[j])[1]
                

        chromosome_matrix_prev <- chromosome_matrix

}

best_chromosome <- chromosome_matrix[,elitist_i[1]]
best_index <- which(best_chromosome == 1)
best_index <- append(best_index,ncol(train))

model <- lm(b.p. ~ ., data= train[,best_index])
print(summary(model))
plot(F,type="l", xlab="Generation") 

################Print plots
train <- train[,best_index]
test <- test[,best_index]

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
test_R2 <- 1 - (1-test_R2)*(nrow(test)-1)/(nrow(test) - (ncol(test)-1) -1)                                                                                                                     
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
