## Case Study 3
## TODO Smooth mortality data used in model fitting, forecasting and simulations
## You can find further instructions commented out in the code
library(demography)
library(StMoMo)
library(rgl)
library(googleVis)
library(fanplot)

keep(skDemo, BEL, SCR, BEL_const, BEL_case1, SCR_case1, BEL_case2, SCR_case2, df_2gether, sure = TRUE)

source("temp_functions.R")

load("skDemo.RData")

## TODO using the help smooth the mortality data by adjusting the commented code
# ?smooth.demogdata
# attributes(skDemo)
# skDemo <- 
# years <- 
# ages<- 
# Dxt <- 
# E0xt <- 
# Ecxt <- 
# qxt <- 

forecastTime <- 120
ages.fit <- 25:90

persp3d(ages[0:100], years, qxt[0:100,], col="skyblue", shade=TRUE,xlab="Ages (0-100)",
        ylab="Years",zlab="Mortality rate (log)")

weights <- genWeightMat(ages.fit, years, 3)

## Modelling

modelFit <- array(data = NA, c(7, 2))
colnames(modelFit) <- c("AIC", "BIC")
rownames(modelFit) <- c("LC", "APC", "CBD", "M6", "M7", "M8", "PLAT")

## LC model under a Binomial setting - M1
LC <- lc(link = "logit")
LCfit <- fit(LC, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

plot(LCfit, nCol = 3)

LCres <- residuals(LCfit)
plot(LCres, type = "scatter", reslim = c(-3.5, 3.5))

modelFit[1, 1] <- AIC(LCfit)
modelFit[1, 2] <- BIC(LCfit)

# Forecast
LCfor <- forecast(LCfit, h = forecastTime)
LCqxt <- cbind(LCfor$fitted, LCfor$rates)

## APC Currie (2006) - M3
APC <- apc("logit")
APCfit <- fit(APC, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights, start.ax = LCfit$ax,
              start.bx = LCfit$bx, start.kt = LCfit$kx)

modelFit[2, 1] <- AIC(APCfit)
modelFit[2, 2] <- BIC(APCfit)

APCfor <- forecast(APCfit, h = forecastTime, gc.order = c(1, 1, 0))
APCqxt <- cbind(APCfor$fitted, APCfor$rates)

## CBD model Cairns (2009) under a Binomial distribution of deaths Haberman and Renshaw (2011) - M5
CBD <- cbd()
CBDfit <- fit(CBD, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[3, 1] <- AIC(CBDfit)
modelFit[3, 2] <- BIC(CBDfit)

CBDfor <- forecast(CBDfit, h = forecastTime)
CBDqxt <- cbind(CBDfor$fitted, CBDfor$rates)

## reduced PLAT M6
M6 <- m6()
M6fit <- fit(M6, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)
modelFit[4, 1] <- AIC(M6fit)
modelFit[4, 2] <- BIC(M6fit)

M6for <- forecast(M6fit, h = forecastTime, gc.order = c(2, 0, 0))
M6qxt <- cbind(M6for$fitted, M6for$rates)

## M7 under Binomial setting
M7 <- m7(link = "logit")
M7fit <- fit(M7, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[5, 1] <- AIC(M7fit)
modelFit[5, 2] <- BIC(M7fit)

M7for <- forecast(M7fit, h = forecastTime, gc.order = c(2, 0, 0))
M7qxt <- cbind(M7for$fitted, M7for$rates)

## M8
M8 <- m8(link = "logit", xc = 85)
M8fit <- fit(M8, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)
modelFit[6, 1] <- AIC(M8fit)
modelFit[6, 2] <- BIC(M8fit)

M8for <- forecast(M8fit, h = forecastTime, gc.order = c(2, 0, 0))
M8qxt <- cbind(M8for$fitted, M8for$rates)

## PLAT
f2 <- function(x, ages) mean(ages) - x
constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages){
  nYears <- dim(wxt)[2]
  x <- ages
  t <- 1:nYears
  c <- (1 - tail(ages, 1)):(nYears - ages[1])
  xbar <- mean(x)
  #\sum g(c)=0, \sum cg(c)=0, \sum c^2g(c)=0
  phiReg <- lm(gc ~ 1 + c + I(c^2), na.action = na.omit)
  
  phi <- coef(phiReg)
  gc <- gc - phi[1] - phi[2] * c - phi[3] * c^2
  kt[2, ] <- kt[2, ] + 2 * phi[3] * t
  kt[1, ] <- kt[1, ] + phi[2] * t + phi[3] * (t^2 - 2 * xbar * t)
  ax <- ax + phi[1] - phi[2] * x + phi[3] * x^2
  #\sum kt[i, ] = 0
  ci <- rowMeans(kt, na.rm = TRUE)
  ax <- ax + ci[1] + ci[2] * (xbar - x)
  kt[1, ] <- kt[1, ] - ci[1]
  kt[2, ] <- kt[2, ] - ci[2]
  list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
}
PLAT <- StMoMo(link = "logit", staticAgeFun = TRUE,
               periodAgeFun = c("1", f2), cohortAgeFun = "1",
               constFun = constPlat)
PLATfit <- fit(PLAT, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[7, 1] <- AIC(PLATfit)
modelFit[7, 2] <- BIC(PLATfit)

PLATfor <- forecast(PLATfit, h = forecastTime, gc.order = c(2, 0, 0))
PLATqxt <- cbind(PLATfor$fitted, PLATfor$rates)

## Collect fitted models for simulation purposes
modelsFitted <- list(LC = LCfit, APC = APCfit, 
                     CBD = CBDfit,  M6 = M6fit, M7 = M7fit, M8 = M8fit, PLAT = PLATfit)

## Model fit criteria AIC and BIC
modelFit

## Ploting - inspection of forecasted mortality rates for 65 year old
years_chart <- c(years, (years[length(years)]+1):(years[length(years)]+forecastTime))
plot(years,(Dxt/E0xt)["65", ], pch=21, col='blue', bg='lightblue',
     xlim = range(years_chart), ylim = range(0.01, 0.03), xlab = "Years projection", ylab = "Mortality rate at age 65",
     main = "Forecasts with different models for qxt", bty="n")
lines(years_chart, LCqxt["65",], col = "red", lwd = 2)
lines(years_chart, CBDqxt["65",], col = "lightblue", lwd = 2)
lines(years_chart, APCqxt["65",], col = "green", lwd = 2)
lines(years_chart, M6qxt["65",], col = "pink", lwd = 2)
lines(years_chart, M7qxt["65",], col = "yellow", lwd = 2)
lines(years_chart, M8qxt["65",], col = "purple", lwd = 2)
lines(years_chart, PLATqxt["65",], col = "grey", lwd = 2)
abline(h = 0, v = 2009, col = "gray60")

legend("topright", c("LC", "CBD", "APC", "M6", "M7", "M8", "PLAT"), col=c("red", "lightblue", "green",
                                                                          "pink", "yellow", "purple", "grey"),
       lty=c(1,1), lwd =2, cex = 1, x.intersp=0.7, y.intersp = 1, bty = "n")


## Extrapolate data to omega age = 120
extrapolate <- kannistoExtrapolation(PLATqxt, ages.fit, years_chart)
PLATqxtExtr <- extrapolate$qxt

LCextrapolate <- kannistoExtrapolation(LCqxt, ages.fit, years_chart)
LCqxtExtr <- LCextrapolate$qxt

APCextrapolate <- kannistoExtrapolation(APCqxt, ages.fit, years_chart)
APCqxtExtr <- APCextrapolate$qxt

CBDextrapolate <- kannistoExtrapolation(CBDqxt, ages.fit, years_chart)
CBDqxtExtr <- CBDextrapolate$qxt

M6extrapolate <- kannistoExtrapolation(M6qxt, ages.fit, years_chart)
M6qxtExtr <- M6extrapolate$qxt

M7extrapolate <- kannistoExtrapolation(M7qxt, ages.fit, years_chart)
M7qxtExtr <- M7extrapolate$qxt

M8extrapolate <- kannistoExtrapolation(M8qxt, ages.fit, years_chart)
M8qxtExtr <- M8extrapolate$qxt

models <- list(LCqxtExtr = LCqxtExtr, APCqxtExtr = APCqxtExtr, 
               CBDqxtExtr = CBDqxtExtr,  M6qxtExtr = M6qxtExtr, 
               M7qxtExtr = M7qxtExtr, M8qxtExtr = M8qxtExtr, PLATqxtExtr = PLATqxtExtr)

## Annuity projection
## Assumption annuity due deffered at 65

## Read in portfolio data
portfolio <- read.csv('portfolio.csv')

## Read in experience data
experience.factors <- read.csv('experience-factors.csv')

ages.fit <- 25:120
valyear <- 2015
pension <- 1000

# calculate ages of insured
portfolio$age <- valyear - portfolio$YoB
experience.factors$total <- (experience.factors$Male + experience.factors$Female)/2

expF <- experience.factors$total[ages.fit]

BEL_case3 <- array(NA, c(7,1))
rownames(BEL_case3) <-  c("LC", "APC", "CBD", "M6", "M7", "M8", "PLAT")
colnames(BEL_case3) <- "BEL"


for (m in 1:length(models)){
  output <- list()
  output2 <- list()
  for (i in 1:nrow(portfolio)){
    output[[i]] <- DFcashflow(models[[m]]*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, valyear = valyear, ir = 0.02, type = 1)*pension
    output2[[i]] <- DFcashflow(models[[m]]*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, valyear = valyear, ir = 0.02, type = 2)*portfolio$Premium[i]
  }
  BEL_case3[m, 1] <- do.call(sum, output)-do.call(sum, output2)
}


## Show BEL per model
df_case3 = data.frame(models =  c("LC", "APC", "CBD", "M6", "M7", "M8", "PLAT"), BEL = BEL_case3)
Column_case3 <- gvisColumnChart(df_case3, options=list(series="[{color:'#0A79BF'}]"))
plot(Column_case3)

## SIMULATIONS for SCR calculation purposes
set.seed(1234)

# number of simulations
nsim <- 200
models2run <- length(modelsFitted)
ages.fit <- 25:90

modelSim <- list()
selectBEL <- list()

for (m in 1:7){
  modelSim[[m]] <- simulate(modelsFitted[[m]], nsim = nsim, h = forecastTime)
  collectBEL <- array(NA, c(200,1))
  for (s in 1:nsim){
    prem_s <- 0
    ben_s <- 0
    qx <- cbind(modelSim[[m]]$fitted[, , s], modelSim[[m]]$rates[, , s])
    extrapolate <- kannistoExtrapolation(qx, ages.fit, years_chart)
    for (i in 1:nrow(portfolio)){
      ben_s <- ben_s + DFcashflow(extrapolate$qxt*expF, ageStart = portfolio$age[i], omegaAge = 120, 
                                  pensionAge = 65, valyear = valyear, ir = 0.02, type = 1)*pension
      prem_s <- prem_s + DFcashflow(extrapolate$qxt*expF, ageStart = portfolio$age[i], omegaAge = 120, 
                                    pensionAge = 65, valyear = valyear, ir = 0.02, type = 2)*portfolio$Premium[i]
    }
    collectBEL[s, 1] <- ben_s - prem_s
  }
  selectBEL[[m]] <- quantile(collectBEL, probs = 0.995, type = 1)
}

SCR_case3 <- as.numeric(selectBEL) - BEL_case3
colnames(SCR_case3) <- "SCR"

## Plot simulations for LC model
qxt <- Dxt/E0xt
plot(LCfit$years, qxt["65", ], xlim = c(1950, 2129), ylim = range(PLATqxt["65",]),
     xlab = "Years", ylab = "Mortality rates", main = "Mortality rates at age 65",
     pch = 20, log = "y", type = "l")
matlines(modelSim[[1]]$years, modelSim[[1]]$rates["65", , 1:20], type = "l", lty = 1, col = 1:20)

## Plot model uncertainity
probs <- c(0.5, 2.5, 10, 25, 75, 90, 97.5, 99.5)
plot(LCfit$years, qxt["65", ], xlim = c(1950, 2129), ylim = c(0.0025, 0.04),
     xlab = "Years", ylab = "Mortality rate", main = "Uncertainity associated with a model forecast",
     pch = 20, log = "y")
fan(t(modelSim[[1]]$rates["65", , ]), start = 2010, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("yellow", "darkgreen")), ln = NULL)

## Show SCR per model
df_case3 = data.frame(models =  c("LC", "APC", "CBD", "M6", "M7", "M8", "PLAT"), BEL = BEL_case3, SCR = SCR_case3)
Column_case3 <- gvisColumnChart(df_case3,
                          options=list(series="[{color:'#0A79BF'}, {color: '#B2246B'}]"))
plot(Column_case3)

## TODO Plot results from all cases 