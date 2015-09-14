## Case Study 0
## Constant mortality rates
# library(demography)
# library(StMoMo)
# library(rgl)
# library(googleVis)
# library(fanplot)

load("skDemo.RData")
years <- skDemo$year
ages<- skDemo$age
Dxt <- skDemo$rate$total * skDemo$pop$total
E0xt <- skDemo$pop$total + 0.5 * Dxt
Ecxt <- skDemo$pop$total
qxt <- Dxt/E0xt

# TODO create the qxt table for ages from 25 to 90 and years from 1950 to 2129
# Next set mortality rates constant  to mortalities from 2009 for years 2010 to 2129
head(qxt)
dim(qxt)
?Extract
ages.fit <- 25:90
CONSTqxt <- qxt[ages.fit, pmin(1:(ncol(qxt)+120), ncol(qxt))]
rownames(CONSTqxt) <- as.character(ages.fit)
colnames(CONSTqxt) <- as.character(1950:2129)

## TODO peform extrapolation to 120
extrapolate <- kannistoExtrapolation(qx = CONSTqxt, ages = ages.fit, years = 1950:2129)
CONSTqxtExtr <- extrapolate$qxt

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

#Calculate BEL for the portfolio
#TODO adjust the code below
output <- list()
output2 <- list()
for (i in 1:nrow(portfolio)){
  output[[i]] <- DFcashflow(CONSTqxtExtr*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, valyear = valyear, ir = 0.02, type = 1)*pension
  output2[[i]] <- DFcashflow(CONSTqxtExtr*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, valyear = valyear, ir = 0.02, type = 2)*portfolio$Premium[i]
}
BEL_const <- do.call(sum, output)-do.call(sum, output2)

## Show BEL per model
df = data.frame(models = c("LC", "APC", "CBD", "M6", "M7", "M8", "PLAT", "CONST"), BEL = c(BEL, BEL_const))
Column <- gvisColumnChart(df, options=list(series="[{color:'#0A79BF'}]"))
plot(Column)
