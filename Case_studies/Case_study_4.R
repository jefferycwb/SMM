## Case Study 4 #constant mortality rates
# TODO create the qxt table for ages from 25 to 90 and years from 1950 to 2129
# And set mortality rates constant  to mortalities from 2009 for years 2010 to 2129

## Note that we already have data loaded in our work space

## First let's understand the data you are working with. 
## Let's have a look and execute the following commands <Ctrl+Enter> runs one line, where cursor is placed
# head(qxt)
# dim(qxt)
# head(ages.fit)
# x<-1:100 
# x[1,length(x),1] #first element, last element and last again
## Check the following functions, they are usefull to solve the case study
# ?Extract ## Hint look into Help examples
# ?range

## First in the console try to extract first column of qxt table, then the last one.
## Try to extract more columns at once?
## This should give you the feeling how R works. Now let's jump into solving the case.

# CONSTqxt <- 
# rownames(CONSTqxt) <- as.character(<put the row names here>)
# colnames(CONSTqxt) <- as.character(<put the column names here>)

  
## TODO peform extrapolation to 120
# extrapolate <- kannistoExtrapolation(qx = , ages = , years = )
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

##Calculate BEL for the portfolio
##TODO adjust the code below
output <- list()
output2 <- list()
# for (i in 1:nrow(portfolio)){
#   output[[i]] <- DFcashflow(<put the right table here>*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, valyear = valyear, ir = 0.02, type = 1)*pension
#   output2[[i]] <- DFcashflow(<put the right table here>*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, valyear = valyear, ir = 0.02, type = 2)*portfolio$Premium[i]
# }
BEL_const <- do.call(sum, output)-do.call(sum, output2)

## Show BEL per model
df = data.frame(models = c("LC", "APC", "CBD", "M6", "M7", "M8", "PLAT", "CONST"), BEL = c(BEL, BEL_const))
Column <- gvisColumnChart(df, options=list(series="[{color:'#0A79BF'}]"))
plot(Column)
