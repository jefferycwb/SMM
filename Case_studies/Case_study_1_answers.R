## CASE STUDY 1

library(googleVis)

# #TODO # Adjust the portfolio data and make everybody + 5 years older
portfolio$age <- portfolio$age + 5

ages.fit <- 25:120
valyear <- 2015
pension <- 1000

BEL_case1 <- array(NA, c(7,1))
rownames(BEL_case1) <- c("LC", "APC", "CBD", "M6", "M7", "M8", "PLAT")
colnames(BEL_case1) <- "BEL"

for (m in 1:length(models)){
  output <- list()
  output2 <- list()
  for (i in 1:nrow(portfolio)){
    output[[i]] <- DFcashflow(models[[m]]*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, valyear = valyear, ir = 0.02, type = 1)*pension
    output2[[i]] <- DFcashflow(models[[m]]*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, valyear = valyear, ir = 0.02, type = 2)*portfolio$Premium[i]
  }
  BEL_case1[m, 1] <- do.call(sum, output)-do.call(sum, output2)
}

## Show BEL per model
df_case1 = data.frame(models = c("LC", "APC", "CBD", "M6", "M7", "M8", "PLAT"), BEL = BEL_case1)
Column_case1 <- gvisColumnChart(df_case1, options=list(series="[{color:'#0A79BF'}]"))
plot(Column_case1)

## SIMULATIONS for SCR calculation purposes
set.seed(1234)

# number of simulations
nsim <- 200
models2run <- length(modelsFitted)
ages.fit <- 25:90
selectBEL <- list()

for (m in 1:models2run){
  #modelSim[[m]] <- simulate(modelsFitted[[m]], nsim = nsim, h = forecastTime)
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

SCR_case1 <- as.numeric(selectBEL) - BEL_case1
colnames(SCR_case1) <- "SCR"

## Plot simulations for PLAT model
qxt <- Dxt/E0xt
plot(PLATfit$years, qxt["65", ], xlim = c(1960, 2129), ylim = range(PLATqxt["65",]),
     xlab = "Years", ylab = "Mortality rates", main = "Mortality rates at age 65",
     pch = 20, log = "y", type = "l")
matlines(modelSim[[1]]$years, modelSim[[1]]$rates["65", , 1:20], type = "l", lty = 1, col = 1:20)

## Plot model uncertainity
probs <- c(0.5, 2.5, 10, 25, 75, 90, 97.5, 99.5)
plot(PLATfit$years, qxt["65", ], xlim = c(1960, 2129), ylim = c(0.0025, 0.04),
     xlab = "Years", ylab = "Mortality rate", main = "Uncertainity associated with a model forecast",
     pch = 20, log = "y")
fan(t(modelSim[[1]]$rates["65", , ]), start = 2010, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("yellow", "darkgreen")), ln = NULL)

## Show BEL per model
df_case1 = data.frame(models = c("LC", "APC", "CBD", "M6", "M7", "M8", "PLAT"), BEL = BEL_case1, SCR = SCR_case1)
Column_case1 <- gvisColumnChart(df_case1,
                          options=list(series="[{color:'#0A79BF'}, {color: '#B2246B'}]"))
plot(Column_case1)

# #TODO create a plot showing BEL and SCR from main run vs. Case 1
# Hint use the code above and check the Help for googleVis function
?gvisColumnChart

df_2gether = data.frame(models = c("LC", "APC", "CBD", "M6", "M7", "M8", "PLAT"), BEL_org = BEL, SCR_org = SCR, BEL = BEL_case1, SCR = SCR_case1)
colnames(df_2gether) <- c("Model", "BEL", "SCR", "BEL case 1", "SCR case 1")
Column_2gether <- gvisColumnChart(df_2gether)
plot(Column_2gether)
