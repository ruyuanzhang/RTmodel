function paramsfitmulti=fitRTmodel_optimize(data)
% fit the RT model. The <data> is a matrix with size nTrial x data
% conditions. The first two columns must be accuracy and RT
% 
%
if notDefined('data')
    error('Please input data');
end

%% some high-level settings
nFit = 100; % fit the model how many times and get the maximum likelihood

%% massage data
nTrials = size(data, 1);
nExpVars = size(data, 2)-2; % the first two colume are accuracy and RT
% extract experiment conditions
cond = cell(1,nExpVars);
nCond = 1;
for i = 1:nExpVars
    cond{i} = unique(data(:,2+i));
    nCond = nCond * length(cond{i});
end
% extract data and form a data cell. A data cell is easier to fit the model.
% This is ugly..
dataSorted = cell(1,nCond);
tmpidx=1;
for i = 1:length(cond{1})
    for j=1:length(cond{2})
        for k=1:1:length(cond{3})
            idx = find(data(:,3)==cond{1}(i) & data(:,4)==cond{2}(j) & data(:,5)==cond{3}(k));
            dataSorted{tmpidx}=[data(idx,1) data(idx,2)];
            tmpidx=tmpidx+1;
        end
    end
end
clear tmpidx idx
% dataSorted is a 1 X nCond cell that contains data for all experimental
% conditions
nParams = 2 + 2 * nCond;  % number of free parameters
minRT = min(data(:,2));
maxRT = max(data(:,2));

%% Do the fitting, we used EM algorithm to iteratively approaching the model fitting
% optimize setting
LB = [1e-5, 1e-5, 1 * ones(1, 2 * nCond)];  
UB = [10, minRT-1, maxRT * ones(1, 2* nCond)]; % alpha, phi, lambda_tim (8), lambda_rim (8)
options = optimset('Algorithm','sqp','MaxFunEvals',1e+5,'MaxIter',1e+5);

% do it
paramsfitmulti = zeros(nFit,length(LB)+3); % extra 3 spots should include posloglikelihood, AIC and BIC values
fprintf('In the model fitting process...wait...\n')
for i = 1:nFit % loop fitting
    progressbar(i,nFit);
    params0 = rand(1,nParams).*(UB-LB)+LB; % give a random initial values
    % fit
    [paramsfitmulti(i,1:end-3), paramsfitmulti(i,end-2)]= fminsearchbnd(@(params) computeposloglikeli(params,dataSorted), params0,LB, UB, options);
    % calculate AIC and BIC
    paramsfitmulti(i,end-1) = 2*nParams + 2*paramsfitmulti(i,end-2); % AIC
    paramsfitmulti(i,end) = log(nTrials)*nParams + 2*paramsfitmulti(i,end-2); % BIC
end
fprintf('Model fitting done !\n');

end


function posloglikeli=computeposloglikeli(params, data)
% This is the function for calculating the positive loglikelihood given one
% subject's data. the RT model is derived from the ref below.
%
% Glickman, M. E., Gray, J. R., & Morales, C. J. (2005). 
% Combining speed and accuracy to assess error-free cognitive processes. psychometrika, 70(3), 405-425.
%
% Input:
%   <params>: 1x18 vector. The 1st and 2nd are alpha and phi. The 3rd to
%       10th are eight lambda_tim parameters for eight conditions.
%       Similarly, the 11th to 18th are eight lambda_rim values

%   <data>: a 1 x conditions cell, a data of single analysis
%
% Output:
%   <posloglikeli>: the postive loglikelihood to minize

%% some set up 
nParams = length(params); 
alpha = params(1);
phi = params(2);
lambda_tim = params(3:3+(nParams-2)/2-1);
lambda_rim = params(3+(nParams-2)/2:end);
chance = 0.5;% note that this depends on the task, chance level for a 2AFC task is 0.5

%% some simple calculation
t_lambda_tim = zeros(size(lambda_tim)); % intermediate paramters according to the paper
t_lambda_rim = zeros(size(lambda_rim));

posloglikeli = 0;
for iCond = 1:numel(data)
    
    data_cond = data{iCond};
    
    nTrials = size(data_cond,1);
    correct = data_cond(:,1);
    RT = data_cond(:,2);
    n1 = sum(correct); % number of correct trials
    n0 = nTrials-n1; % number of incorrect trials
    
    
    % convert to parameters described in the paper. This is important since
    % it is much easier to set the parameter range for optimization
    t_lambda_tim(iCond) = lambda_tim(iCond)^(-alpha);
    t_lambda_rim(iCond) = lambda_rim(iCond)^(-alpha);
    
    %% compute pos loglikeli
    % We calculated the positive loglikelihood. Maximizing loglikelihood is
    % equivalent to minimzing postive loglikelihood
    % This calculation is based the Equation (7) in the reference.
    term1 = nTrials*log(alpha/(t_lambda_tim(iCond) + t_lambda_rim(iCond))) + 2*n0*log(t_lambda_rim(iCond)) + n0*log((1-chance)) + n1*log((t_lambda_tim(iCond)^2+t_lambda_rim(iCond)^2*chance));
    term2 = sum((alpha-1)*log(RT-phi) - (t_lambda_tim(iCond) + t_lambda_rim(iCond))*(RT-phi).^alpha);
    
    posloglikeli = posloglikeli - (term1 + term2); % We add negative to convert negative loglikeli to posloglikeli;
    
end

end