function paramsfitmulti=fitRTmodel(data)
%
%
%
if notDefined(data)
    error('Please input data');
end

%% some high-level settings
maxIter = 10000; % maximum iteration
chance = 0.5; % chance level of this task

%% massage data
nTrials = size(data, 1);
nExpVars = size(data, 2)-2; % the first two colume are accuracy and RT
% extract experiment conditions
cond = cell(1,nExpVars);
nCond = 1;
for i = 1:nExpVars
    cond{i} = unique(data(:,2+i));
    nCond = nCond * cond{i};
end
% extract data and form a cell, this is ugly..
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
% dataSorted is a 1xnCond cell that contains data for all experiment
% conditions

nParams = 2 * nCond + 2; % number of parameters to estimate
%% Do the fitting, we used EM algorithm to iteratively approaching the model fitting
% optimize setting
LB = 1e-5 * ones(1,nParams);  % lambdat, lambdar, alpha, phi
UB = 100 * ones(1,nParams);
options = optimset('Algorithm','active-set','MaxFunEvals',1e+4,'MaxIter',1e+4);

% do it
fprintf('In the model fitting process...wait...\n')
for i = 1:nFit % loop fitting
    i
    params0 = [UB(1)*rand, rand*(UB(2)-LB(2))+LB(2), UB(3)*rand, UB(4)*rand];
    [paramsfitmulti(i,1:4), paramsfitmulti(i,5)]= fminsearchbnd(@(params) computeposloglikeli(params,data), params0,LB, UB, options);
    paramsfitmulti(i,6) = 2*nParams + 2*paramsfitmulti(i,5); % AIC
    paramsfitmulti(i,7) = log(nTrials)*nParams + 2*paramsfitmulti(i,5); % BIC
end
% derive minimal likeli


%%
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
%   <params>: 1x4 vector, the four params being estimated from the fitting
%       lambdat: 
%       lambdar:
%       alpha:
%       phi:
%   <data>: a 2xtrial matrix, a data of single analysis
%
% Output:
%   <posloglikeli>: We 

%% some set up 
lambdat = params(1);
lambdar = params(2);
alpha = params(3);
phi = params(4);

chance = 0.5;% note that this depends on the task, chance level for a 2AFC task is 0.5

%% some simple calculation
nTrials = size(data,2);
correct = data(1,:);
RT = data(2,:);

n1 = sum(correct); % number of correct trials
n0 = nTrials-n1; % number of incorrect trials

%% compute pos loglikeli
% we calculated the positive loglikelihood. Maximizing loglikelihood is
% equivalent to minimzing postive loglikelihood
% This calculation is based the Equation (7) in the paper.
% 
term1 = (alpha/(lambdat + lambdar))^nTrials*lambdar^(2*n0)*(lambdat^2+lambdar^2*chance)^n1;
term2 = prod((RT-phi)^(alpha-1)*exp(-((lambdat + lambdar)*(RT-phi)^alpha)));

posloglikeli = -term1*term2; % We add negative to convert loglikeli to posloglikeli;



end