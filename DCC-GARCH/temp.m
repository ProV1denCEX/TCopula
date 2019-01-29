cSetupPlatform.Allocate_Setting.AdjStep = 0;
cSetupPlatform.Allocate_Setting.RiskFreeRate = 0;
cSetupPlatform.Allocate_Setting.AdjInterval = 10;
cSetupPlatform.Allocate_Setting.Date_Start = 20171001;
cSetupPlatform.Allocate_Setting.Date_End = 20181028;
cSetupPlatform.Allocate_Setting.SampleDays = 150;

cSetupPlatform.Allocate_Setting.AllocateMethod = {'Equal Weight',...
    'Risk Parity - CVaR - History', ...
    'Risk Parity - CVaR - Norm', ...
    'Risk Parity - CVaR - T'};

cSetupPlatform.Markowitz_Setting.GivenRisk = 0.05;

cSetupPlatform.Simu.Road = 500;
cSetupPlatform.Simu.Len = cSetupPlatform.Allocate_Setting.AdjInterval;
cSetupPlatform.Simu.InterLag = 0.1;


%% Get the Checkpoint
dTimeAxis = cNv2Allocate(1).Data(:, 1);
cSetupPlatform.Allocate_Setting.CheckPoint = find(mod(1 : length(dTimeAxis) , cSetupPlatform.Allocate_Setting.AdjInterval) == 0);

%% Generate Portfolio Weight
for iIndexTime = cSetupPlatform.Allocate_Setting.SampleDays : length(cNv2Allocate(1).Data) - 1
    disp(['Generating Portfolio ', num2str(iIndexTime), ' / ', num2str(length(cNv2Allocate(1).Data) - 1), ' Finished'])
    if ismember(iIndexTime, cSetupPlatform.Allocate_Setting.CheckPoint) ||...
            iIndexTime == cSetupPlatform.Allocate_Setting.SampleDays
        
        [dMean, dCov, dPool, dRev] = Fun_Cal_MeanCov(iIndexTime);
        for iMethod = 1 : length(cAllocateResult)
            sMethod = cSetupPlatform.Allocate_Setting.AllocateMethod{iMethod};
            
            switch sMethod
                case 'Equal Weight'
                    Fun_Allocate_EW(iIndexTime, iMethod, dCov, dPool);

                case 'Risk Parity - CVaR - History'
                    Fun_Allocate_RPC(iIndexTime, iMethod, dCov, dRev, dPool);
                    
                case 'Risk Parity - CVaR - Norm'
                    Fun_Allocate_RPC(iIndexTime, iMethod, dCov, dRev, dPool);
                    
                case 'Risk Parity - CVaR - T'
                    Fun_Allocate_RPC(iIndexTime, iMethod, dCov, dRev, dPool);
                    
                otherwise
            end
        end
    else
        for iMethod = 1 : length(cSetupPlatform.Allocate_Setting.AllocateMethod)
            cAllocateResult(iMethod).Result(iIndexTime + 1, 2 : end) = cAllocateResult(iMethod).Result(iIndexTime, 2 : end);
        end
    end
end


%% ****************************** Function **************************************
function Fun_Allocate_RPC(iIndexTime, nSerial, dCov, dRev, dPool)
global cAllocateResult
global cSetupPlatform

dWeigth_Last = cAllocateResult(nSerial).Result(iIndexTime, 2 : end);

% Update Return
cSetupPlatform.Markowitz_Setting.Return = dRev;


% Choose the Optimize Target
sMethod = cAllocateResult(nSerial).Method;
switch sMethod
    case 'Risk Parity - CVaR - History'
        % 历史分布拟合
        sTarget = @Fun_Target_CVaR_Hist;
        
        % Cov - Achieved
        cSetupPlatform.Markowitz_Setting.Cov = dCov;
        
    case 'Risk Parity - CVaR - Norm'
        % 正态分布拟合
        sTarget = @Fun_Target_CVaR_Norm;
        
        % Cov - DCC(1, 1)
        dTemp = dRev;
        dMean = mean(dTemp);
        dMean = repmat(dMean, size(dTemp, 1), 1);
        dTemp = dTemp - dMean;
        [~, ~, dCov] = dcc(dTemp,[],1,0,1);
        
        cSetupPlatform.Markowitz_Setting.Cov = dCov(:, :, end);
        
    case 'Risk Parity - CVaR - T'
        % T 分布拟合
        sTarget = @Fun_Target_CVaR_T;
        
         % Cov - DCC(1, 1)
        dTemp = dRev;
        dMean = mean(dTemp);
        dMean = repmat(dMean, size(dTemp, 1), 1);
        dTemp = dTemp - dMean;
        [~, ~, dCov] = dcc(dTemp,[],1,0,1);
        
        cSetupPlatform.Markowitz_Setting.Cov = dCov(:, :, end);
        
    otherwise
        
end

% Given Risk, Parity CVaR
dConstraint_A = [];
dConstraint_b = [];
sConstraint_Nonl = @Fun_Constraint_NonL;

% Check the Limit
[dConstraint_X0, dConstraint_vlb, dConstraint_vub, dConstraint_Aeq, dConstraint_beq] = ...
    Fun_Get_Constraint(iIndexTime, dWeigth_Last, dPool);

% Optimization
dWeightAllocate = fmincon( ...
    sTarget, ...
    dConstraint_X0, ...
    dConstraint_A, dConstraint_b,...
    dConstraint_Aeq, dConstraint_beq, ...
    dConstraint_vlb, dConstraint_vub, ...
    sConstraint_Nonl);

% Check the Step
dWeight = Fun_Check_AdjStep(dWeightAllocate);

% Write
cAllocateResult(nSerial).Result(iIndexTime + 1, dPool + 1) = dWeight;

end

function nOutput = Fun_Target_CVaR_Hist(dWeight)
global cSetupPlatform

% Cal CVaR
dReturnFull_Sample = cSetupPlatform.Markowitz_Setting.Return;
dSimuRoot = [(1 : length(dReturnFull_Sample))', sort(dReturnFull_Sample)];
dCVaR = zeros(size(dSimuRoot, 2) - 1, 1);
for iData = 1 : size(dReturnFull_Sample, 2)
    dTemp = Fun_Get_Interpolation(dSimuRoot(:, [1, iData + 1]), cSetupPlatform.Simu.InterLag, 'Spline');
    dLocated = dTemp(:, 2) <= prctile(dTemp(:, 2), 0.05);
    dCVaR(iData) = mean(dTemp(dLocated, 2));
end

% 2 Parity
dCVaR_All = dCVaR .* dWeight(:);
nPortRisk = sum(dCVaR_All);
nOutput = sum((dCVaR_All / nPortRisk - nPortRisk / length(dWeight(:))).^2);
end

function nOutput = Fun_Target_CVaR_Norm(dWeight)
global cSetupPlatform

% Cal CVaR
dReturnFull_Sample = cSetupPlatform.Markowitz_Setting.Return;
dSimuRoot = [(1 : length(dReturnFull_Sample))', sort(dReturnFull_Sample)];
dCVaR = zeros(size(dSimuRoot, 2) - 1, 1);

% Return Simu
dMean = mean(dReturnFull_Sample);
dCov = cov(dReturnFull_Sample);
dReturnSimu = mvnrnd(dMean, dCov, cSetupPlatform.Simu.Len * cSetupPlatform.Simu.Road);

for iData = 1 : size(dReturnFull_Sample, 2)
    dTemp = dReturnSimu(iData);
    % Cal CVaR
    dLocated = dTemp <= prctile(dTemp, 0.05);
    dCVaR(iData) = mean(dTemp(dLocated));
end

% 2 Parity
dCVaR_All = dCVaR .* dWeight(:);
nPortRisk = sum(dCVaR_All);
nOutput = sum((dCVaR_All / nPortRisk - nPortRisk / length(dWeight(:))).^2);
end

function nOutput = Fun_Target_CVaR_T(dWeight)
global cSetupPlatform

% Cal CVaR
dReturnFull_Sample = cSetupPlatform.Markowitz_Setting.Return;
dSimuRoot = [(1 : length(dReturnFull_Sample))', sort(dReturnFull_Sample)];
dCVaR = zeros(size(dSimuRoot, 2) - 1, 1);

% Scale the Rev 2 fit
dScale_Temp = zeros(size(dReturnFull_Sample));
for iData = 1 : size(dReturnFull_Sample, 2)
    dScale_Temp(:, iData) = ksdensity(dReturnFull_Sample(:, iData),dReturnFull_Sample(:, iData),'function','cdf');
end
[~, nDF] = copulafit('t', dScale_Temp);

dCov = cov(dReturnFull_Sample);

% Simu
dReturnSimu = mvtrnd(dCov, nDF, cSetupPlatform.Simu.Len * cSetupPlatform.Simu.Road);

% Scale back
for iData = 1 : size(dReturnSimu, 2)
    dReturnSimu(:, iData) = ksdensity(dReturnFull_Sample(:, iData),dReturnSimu(:, iData),'function','icdf');
end

% Cal CVaR
for iData = 1 : size(dReturnFull_Sample, 2)
    dTemp = dReturnSimu(iData);
    dLocated = dTemp <= prctile(dTemp, 0.05);
    dCVaR(iData) = mean(dTemp(dLocated));
end

% 2 Parity
dCVaR_All = dCVaR .* dWeight(:);
nPortRisk = sum(dCVaR_All);
nOutput = sum((dCVaR_All / nPortRisk - nPortRisk / length(dWeight(:))).^2);
end