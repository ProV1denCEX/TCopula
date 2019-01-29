library(data.table)
library(xts)
freturn_df = fread('CTA_Species_Return.csv')
freturn_df$Date = as.Date(freturn_df$Date)
freturn_ts = as.xts(freturn_df)

library(copula)

#["2011-08-04/2016-01-31"]
#["2016-02-01/2018-02-14"]
#["2018-02-15/2018-10-26"]

# before_market_crash
gumbel = xvCopula(gumbelCopula(dim = 11), as.matrix(freturn_ts))
frank = xvCopula(frankCopula(dim = 11), as.matrix(freturn_ts))
joe = xvCopula(joeCopula(dim = 11), as.matrix(freturn_ts))
clayton = xvCopula(claytonCopula(dim = 11), as.matrix(freturn_ts))
normal = xvCopula(normalCopula(dim = 11), as.matrix(freturn_ts))
tvalue = xvCopula(tCopula(dim = 11), as.matrix(freturn_ts))

