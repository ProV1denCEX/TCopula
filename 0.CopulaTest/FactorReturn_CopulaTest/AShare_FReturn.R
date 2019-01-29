library(data.table)
library(xts)
freturn_df = fread('10F_1Vol_all_freturn.csv')
freturn_df$Date = as.Date(freturn_df$Date)
freturn_ts = as.xts(freturn_df)

library(copula)

before_market_crash = freturn_ts["2011-08-04/2016-01-31"]
before_VMom_crash = freturn_ts["2016-02-01/2018-02-14"]
after_YTD = freturn_ts["2018-02-15/2018-10-26"]

# before_market_crash
xvCopula(gumbelCopula(dim = 10), as.matrix(before_market_crash)) #dim = len(factor_list)
xvCopula(frankCopula(dim = 10), as.matrix(before_market_crash))
xvCopula(joeCopula(dim = 10), as.matrix(before_market_crash))
xvCopula(claytonCopula(dim = 10), as.matrix(before_market_crash))
xvCopula(normalCopula(dim = 10), as.matrix(before_market_crash))
xvCopula(tCopula(dim = 10), as.matrix(before_market_crash))

# before_VMom_crash
xvCopula(gumbelCopula(dim = 10), as.matrix(before_VMom_crash))
xvCopula(frankCopula(dim = 10), as.matrix(before_VMom_crash))
xvCopula(joeCopula(dim = 10), as.matrix(before_VMom_crash))
xvCopula(claytonCopula(dim = 10), as.matrix(before_VMom_crash))
xvCopula(normalCopula(dim = 10), as.matrix(before_VMom_crash))
xvCopula(tCopula(dim = 10), as.matrix(before_VMom_crash))

# after_YTD
xvCopula(gumbelCopula(dim = 10), as.matrix(after_YTD))
xvCopula(frankCopula(dim = 10), as.matrix(after_YTD))
xvCopula(joeCopula(dim = 10), as.matrix(after_YTD))
xvCopula(claytonCopula(dim = 10), as.matrix(after_YTD))
xvCopula(normalCopula(dim = 10), as.matrix(after_YTD))
xvCopula(tCopula(dim = 10), as.matrix(after_YTD))
