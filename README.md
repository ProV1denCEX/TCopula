# TCopula
Copula Application in Global Market

Work Description:
1.	Copula test -> Why Gaussian copula performed bad
We will investigate the dynamic and asymmetric dependence structure within Chinese equity factors, CTA categories and global markets.
1). For Asian equity market: 10 factor portfolios of Chinese equity market
2). For China commodity market: 11 commodity indices
3). For global equity market: main equity indices
4). For global commodity market: Global CTA Component Index (Bloomberg Commodity Outlook)
For data described above, we shall test Gaussian, student-t, gumbel, frank, joe, clayton copula. We will analyze the out of sample performance given different copulas via leave one out copula information criteria (CIC).

2.	We will choose the copula with largest CIC based on tests in part 1. 
Based on realized tests we have, a copula with fat tails and asymmetric dependence might be better choice rather than other commonly used copulas.

3.	Use Copula-DCC-GARCH to predict OOS joint distribution across different asset classes.
We consider a dynamic asymmetric copula and forecast it for one step joint distribution. For instance, in equity factor copula test, we construct decile factor portfolios sorted on single characteristic: Value(BP_LF), Earnings Yield(1/PE_FTM), Momentum(12-1MOM), 1M Reversal, Volatility(iVol), Size, Beta, Liquidity. Then we predict the joint distribution of these factors’ return by Copula-DCC-GARCH. 
     
4.	Prediction evaluation
There are 3 aspects to measure the performance of predicted copula and Gaussian copula:
1). VaR and CVaR accuracy: 
At every rolling month end, we will generate some scenarios to compute VaR and CVaR of next month. The prediction is based on DCC model fitted on rolling data. Then compare the realized VaR and CVaR next month with our previous prediction. The prediction will be separately computed on Gaussian copula and optimal copula we get in part 2 and 3. Also, the VaR and CVaR will have a volatility regime adjustment to have a better prediction next month.
2). Portfolio performance: 
We will construct a monthly rebalanced portfolio that maximizes CVaR adjusted Alpha returns under sector-neutral and long only constraints. Then we compare rolling OOS portfolio performance with benchmark.
3). Bias Test
We use distance-based statistics to measure predictions bias on correlation matrix across different copulas.
