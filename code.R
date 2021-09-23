### Argus coding Challenge
#0) load the packages.
# load libraries
library(gamlss)
library(gamlss.add)
library(gamlss.dist)
library(corrplot)
library(tseries)
library(DT)
library(dplyr)
library(stats)
library(psych)
library(ggpubr)
library(roll)
library(magrittr)
library(labelVector)


#1)	Take as raw inputs to the data preparation process, the oil data from the gamlss package.
data_oil <- gamlss.data::oil
plot(OILPRICE~SPX_log, data=oil)
#2) Develop a process that allows us to add additional drivers which are transformations of the raw input timeseries. Include the following transformations:
#a. Rolling standard deviation (of arbitrary window)
#b. Rolling mean (of arbitrary window)
#c. Lagging (of arbitrary order)
#d. Leading (of arbitrary order)
#e. Differencing
#f. Spread (between two input drivers)
#g. Ratio (between two input drivers)
#h. Product (between two input drivers)


#window size is 5
data <- as.data.frame(data_oil$OILPRICE)
# assign input to data_oil_2 as a matrix
data_oil_2 <- as.matrix(data_oil$OILPRICE)
# calculate rtolling standard deviation, window = 5
roll_std_deviation <- roll::roll_sd(data_oil_2, 5)
data$roll_std_deviation <- roll_std_deviation#(Rolling standard deviation)
# Rolling mean  
roll_mean <- roll::roll_mean(data_oil_2, 5)
data$roll_mean <- roll_mean# (Rolling mean)
# Lagging  order = 1
data$lag_1 <- dplyr::lag(data_oil$OILPRICE)
# Leading  order = 1
data$lead <- dplyr::lead(data_oil$OILPRICE)
# Differencing
Diff <- data_oil$OILPRICE %>% diff()
Diff[1000] <- NA
data$diff <- Diff
data_oil_transformed <- data
head(data_oil_transformed, n =10)

#f. Spread (between two input drivers) : I am not sure how to calculate Spread between 2 variables. Need little brief about this function.
input=data_oil_transformed[, c("roll_std_deviation", "roll_mean")]
data_temporary <- as.data.frame(input)
data_temporary$Ratio <- data_temporary[,1]/data_temporary[,2]
data_temporary$Product <- data_temporary[,1] * data_temporary[,2]

#3) We must be able to have composition of transformations. Example: First calculate the difference between OILPRICE and resp_LAG, and then calculate the rolling standard deviation.
input=data_oil[,c("OILPRICE", "respLAG")]
data_helper <- as.data.frame(input)
difference <- (data_helper$OILPRICE - data_helper$respLAG)
data_helper$difference <- difference
difference <- as.matrix(difference)
# calculate rtolling standard deviation, window = 7
roll_std <- roll::roll_sd(difference, 7)
data_helper$comp_trans <- roll_std


#4) The sequence of transformations, and which drivers they act on must be specified by the user. One of the main purposes of this challenge is to develop a generic framework to allow this.
#5) For all drivers, either in their raw form or those that results from the application of one or several transformations, we must keep a meta data object where the sequence of transformations is stored. This will allow us to keep track of the meaning of each new driver.

FinalDrivers <- cbind(data_oil_transformed, data_temporary, data_helper)
FinalDrivers <- 
  FinalDrivers[, c("raw_data", "roll_std_deviation", "roll_mean",
                    "lag_1", "lead", "diff", "Ratio", "Product",
                    "comp_trans")]
print_with_label <- function(datarame){
  stopifnot(inherits(datarame, "data.frame"))
  labs <- labelVector::get_label(datarame, names(datarame))
  labs <- sprintf("%s: %s", names(datarame), labs)
  cat("\n")
  cat(labs, sep = "\n")
}
FinalDrivers <-set_label(FinalDrivers,
                          raw_data = "Target",
                          roll_std_deviation = "Rolling standard deviation",
                          roll_mean = "Rolling mean ",
                          lag_1 = "Lagging",
                          lead = "Leading",
                          diff = "Differencing",
                          Ratio = "Ratio",
                          Product = "Multiplication"
)
print_with_label(FinalDrivers)

#6) For each driver that results from the user-specified sequence of transformations, we need to assess a few statistics:
#Normality test
#Stationarity test
#Correlation coefficient with the target
#These statistics need to be stored in the meta data object. The purpose of this is, we may be interested in keeping in the final model only drivers that are normally distributed, or only drivers whose correlation with the target is above a given threshold, or another combination of such criteria.

#a. Normality test Shapiro-Wilk normality test
normality <- function(input_driver) {
  
  print(ggdensity(input_driver, 
                  main = "Density plot of Rolling Standard deviation",
                  xlab = ""))
  
  shapiro.test(input_driver)
  
}

normality(FinalDrivers$roll_std_deviation)

# p-value < 0.05 so the distribution of the data are significantly different from normal distribution. 
#b. Stationarity test Augmented Dickey-Fuller (Adata) t-statistic test hen data is stationary
input=FinalDrivers$roll_mean
input[is.na(input)] <- 0
tseries::adata.test(input)

#c. Correlation coefficient with the target
input=FinalDrivers
input[is.na(input)] <- 0

corrplot(cor(input),
         method = "number",
         type = "upper" # show only upper side
)


##############################################################
#The oil data: fitting GAMLSS models using the RS, CG and mixed algorithms.
#(a) Fit the following model using the RS algorithm 
# 
# m1.rs=gamlss(OILPRICE ~ pb(respLAG) + pb(HO1_log),
#              sigma.formula = ~pb(respLAG) + pb(HO1_log),
#              nu.formula = ~pb(respLAG) + pb(HO1_log),
#              tau.formula = ~pb(respLAG) + pb(HO1_log),
#              family = SHASHo, data = oil, method=RS(20))
# #Since the RS algorithm works, try to fit the model using the mixed algorithm.
# m1.mx=gamlss(OILPRICE ~ pb(respLAG) + pb(HO1_log),
#              sigma.formula = ~pb(respLAG) + pb(HO1_log),
#              nu.formula = ~pb(respLAG) + pb(HO1_log),
#              tau.formula = ~pb(respLAG) + pb(HO1_log),
#              family = SHASHo, data = oil, method=mixed(10,10),
#              gd.tol=Inf)
# 
# GAIC(m1.rs,m1.mx, k=2)
# GAIC(m1.rs,m1.mx, k=3)
# GAIC(m1.rs,m1.mx, k=4)
# GAIC(m1.rs,m1.mx, k=log(nrow(oil)))

#################### ARIMA the to forecast the oil price
ar(input$`data_oil$OILPRICE`)
model = arima(input$`data_oil$OILPRICE`, order = c(20, 3, 1))
nfcs = 20
fcst = predict(model, n.ahead = nfcs)

plot(c(min(input$`data_oil$OILPRICE`), max(input$`data_oil$OILPRICE`) + 
         +            0.3), type = "n", ylab = "SST [C]", xlab = "Date (MM/YY)")
grid()
plot(input$`data_oil$OILPRICE`)
lines(input$`data_oil$OILPRICE`)
pm = fcst$pred
lines(pm, lwd = 2, col = "red")
print('The predicted values are:')
print(pm)
### Due to time limitation and my tough schedule I couldn't tune or try advanced methods. After having the data ready, it is 
#pretty easy to run different methods as they only require some lines of code.