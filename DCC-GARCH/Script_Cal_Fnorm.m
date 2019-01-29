dReturnFull_Sample = cell2mat(cellfun(@(x) {x(1 : end, 2)}, {cNv2Allocate.Data}));
dDate = cNv2Allocate.Data;
dDate = dDate(:, 1);
dNorm = [];
for iIndexTime = 250 : length(dReturnFull_Sample)
    % Get Data
    dReturn_Temp = dReturnFull_Sample(iIndexTime - 249 : iIndexTime, :);
    
    % Half Life
    dRev = diff(dReturn_Temp) ./ dReturn_Temp(1 : end - 1, :);
    dLocated = dRev ~= 0;
    dPool = find(sum(dLocated) >= length(dRev) * 0.9);
    dRev = dRev(:, dPool);
    
    dLocation = floor((1 : cSetupPlatform.RiskParity_Setting.HalfTimes) ...
        * (cSetupPlatform.Allocate_Setting.SampleDays / cSetupPlatform.RiskParity_Setting.HalfTimes)) - 1;
    dRecession = cSetupPlatform.RiskParity_Setting.Recession .^ ...
        ((1 : cSetupPlatform.RiskParity_Setting.HalfTimes) - 1);
    dRecession = dRecession ./ sum(dRecession);
    
    dHCov = zeros(length(dPool), length(dPool));
    for iPart = 1 : cSetupPlatform.RiskParity_Setting.HalfTimes
        if iPart ~= 1
            dRev_Temp = dRev(dLocation(iPart - 1) : dLocation(iPart), :);
        else
            dRev_Temp = dRev(1 : dLocation(iPart), :);
        end
        dCov_Temp = cov(dRev_Temp);
        dHCov = dHCov + dCov_Temp .* dRecession(iPart);
    end
    
    % DCC-Garch-Norm
    dMean = mean(dRev);
    dMean = repmat(dMean, length(dRev), 1);
    dRev_Garch = dRev - dMean;
    [~, ~, dCov_Garch] = dcc(dRev_Garch,[],1,0,1);
    dCov_Garch = dCov_Garch(:, :, end);
    
    % Scale the Rev 2 fit
    dScale_Temp = zeros(size(dRev));
    for iData = 1 : size(dRev, 2)
        dScale_Temp(:, iData) = ksdensity(dRev(:, iData),dRev(:, iData),'function','cdf');
    end
    
    % fit T coupula
    [dCov_Copula, nDF] = copulafit('t', dScale_Temp, 'Method','ApproximateML');
    
    % DCC-Garch-T

    [~, ~, dCov_T] = dcc(dRev_Garch,[],1,0,1, [], [], [], [], [], [], [], [], 'STUDENTST');
    dCov_T = dCov_T(:, :, end);
    
    % Distance
    dNorm(iIndexTime - 249, 1) = norm(cov(dRev), 'fro');
    dNorm(iIndexTime - 249, 2) = norm(dHCov,'fro');
    dNorm(iIndexTime - 249, 3) = norm(dCov_Garch, 'fro');
    dNorm(iIndexTime - 249, 4) = norm(dCov_Copula, 'fro');
    dNorm(iIndexTime - 249, 5) = norm(dCov_T, 'fro');
end

dOutput = [dDate(250 : end), dNorm];
open dOutput