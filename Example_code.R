### Example code to generate point-wise, quantile-based, cluster bootstrap-based, 95% CIs for a GEE model of infection ~ age

#Generate fake dataset
set.seed(100)
practice <- as.data.frame(matrix(NA, nrow=3000, ncol=6))
names(practice) <- c("id", "homeid_cohort", "age", "everinf", "female", "other_var")
practice$id <- seq(1:3000)
practice$homeid_cohort <- round(runif(3000, 1, 300), 0)
practice$age <- round(runif(3000, 2, 14), 0)
practice$everinf <- round(runif(3000, 0, 1), 0)
practice$female <- round(runif(3000, 0, 1), 0)
practice$other_var <- runif(3000, 1, 5000)

#Getting code for legend boxes
source("http://www.math.mcmaster.ca/bolker/R/misc/legendx.R")

#Setting up the color function
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#The model of interest
gee_no2 <- geeglm(everinf ~ age, family=binomial, data=practice, id=homeid_cohort, corstr = "exchangeable", std.err="san.se", scale.fix = TRUE)
summary(gee_no2)


#CI for GEE_no univariate curve

#Preliminaries
practice2 <- dplyr::select(practice, age, id, everinf, homeid_cohort) #subset DF to minimal variables for procedure
practice2 <- practice2[order(practice2$homeid_cohort),] #order the DF by the clustering variable
clusters332 <- names(table(practice2$homeid_cohort)) #vector of names of clustering variable
set.seed(5000) #set seed
full_set <- as.data.frame(matrix(NA, nrow=length(2:14), ncol=10000)) #create DF to hold answers
rownames(full_set) <- 2:14 #age ranges from 2-14 in this case
colnames(full_set) <- paste("Model", seq(1:10000), sep=" ") #model results from 10,000 models
nbs <- 100 #set number of times to do the cluster bootstrap

#Running the cluster bootstrap
for(i in 1:nbs){
  index332 <- sample(1:length(clusters332), length(clusters332), replace=TRUE)
  aa32 <- clusters332[index332]
  bb32 <- table(aa32)
  bootdat32 <- NULL
  for(j in 1:max(bb32)){
    cc32 <- practice2[practice2$homeid_cohort %in% names(bb32[bb32 %in% j]),]
    for(k in 1:j){
      bootdat32 <- rbind(bootdat32, cc32) #bootdat32 is a cluster bootstrap iteration of the dataset, having sampled by the clustering variable
    }
  }
  model_1 <- geeglm(everinf ~ age, family=binomial, data=bootdat32, id=homeid_cohort, corstr = "exchangeable", std.err="san.se", scale.fix = TRUE) #run the model of interest on the iterated dataset
  age_vec <- bootdat32$age #make a vector of the age variable
  prob_vec <- predict(model_1, type="response") #estimate the probability of the outcome based on the model
  answers <- as.data.frame(cbind(age_vec, prob_vec)) #combine the age and prediction vectors
  answers <- unique(answers) #restrict it to unique values of age (ie find the predicted probability for every age based on the model)
  answers <- answers[order(answers$age_vec),] #order the resulting vector by age
  for(l in 1:length(answers$age_vec)){ #transfer the result of this bootstrap iteration into the results vector
    if(answers$age_vec[l] == as.numeric(row.names(full_set)[l])){full_set[l,i] <- answers$V2[l]}
    else(full_set[l,i] <- NULL)
  }
  message(i) #print the run of the bootstrap that was just completed
}
full_set2 <- as.matrix(full_set) #make a matrix for faster processing below
quantiles <- as.data.frame(matrix(NA, nrow=length(2:14), ncol=6)) #make a quantile DF that will be used for plotting
rownames(quantiles) <- 2:14 #name for the rows
names(quantiles) <- c("Lower 90 CI", "Upper 90 CI", "Lower 95 CI", "Upper 95 CI", "Lower 99 CI", "Upper 99 CI") #get 3 kinds of CIs
for(i in 1:length(2:14)){ #find the quantiles of interest
  quantiles[i,1] <- as.numeric(quantile(na.omit(full_set2[i,]), 0.05))
  quantiles[i,2] <- as.numeric(quantile(na.omit(full_set2[i,]), 0.95))
  quantiles[i,3] <- as.numeric(quantile(na.omit(full_set2[i,]), 0.025))
  quantiles[i,4] <- as.numeric(quantile(na.omit(full_set2[i,]), 0.975))
  quantiles[i,5] <- as.numeric(quantile(na.omit(full_set2[i,]), 0.005))
  quantiles[i,6] <- as.numeric(quantile(na.omit(full_set2[i,]), 0.995))
}

#Plot the result
plot(practice$age, practice$everinf, type="n", #Make the plot boundary
     main="Age vs. Probability of infection",
     xlab="Age", ylab="Estimated probability", ylim=c(0, 1), xlim=c(2, 14))
x <- c(2:14, 14:2) #Define the x axis
y90 <- c(quantiles$`Upper 90 CI`, rev(quantiles$`Lower 90 CI`)) #Define the range of the 3 different kinds of CIs
y95 <- c(quantiles$`Upper 95 CI`, rev(quantiles$`Lower 95 CI`))
y99 <- c(quantiles$`Upper 99 CI`, rev(quantiles$`Lower 99 CI`))
salmon90 <- add.alpha("salmon", alpha=0.3) #Define the color scheme for the different kinds of CIs
salmon95 <- add.alpha("salmon", alpha=0.2)
salmon99 <- add.alpha("salmon", alpha=0.1)
polygon(x, y90, col = salmon90, border=NA) #Get a polygon for the CIs across the range of the variable on the x axis
polygon(x, y95, col = salmon95, border=NA)
polygon(x, y99, col = salmon99, border=NA)
lines(practice$age[order(practice$age)], predict(gee_no2, type="response")[order(practice$age)], col="black", lwd=2) #Plot the main curve
salmon90_2 <- add.alpha("salmon", alpha=0.6) #Define the color scheme for the legend's CIs
salmon95_2 <- add.alpha("salmon", alpha=0.3)
salmon99_2 <- add.alpha("salmon", alpha=0.1)
legend(1.3, 1.05, c("Best fitting curve", "90% CI", "95% CI", "99% CI"), bty = "n", x.intersp=0.6, #Add in a legend
       fill=c("black", salmon90_2, salmon95_2, salmon99_2), y.intersp=0.6, box.cex=c(0.6, 0.5))
