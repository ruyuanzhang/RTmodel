% This is the program that fit the RT model. The whole set of program is
% based on the reference below
%

% In this experiment, we have 3 experimental variables and 2^3=6 conditions
% cf: on/off
% cueing: same/diff
% difficulty: hard/easy

% Unlike the paper, we fit the model to individual subject data and
% estimate the parameters. The paper describes four parameters, lamda_tim,
% lambda_rim, phi and alpha. We assume lamda_tim and lambda_rim differ
% across 8 conditions but phi and alpha keep the same. We thus have
% 2*8+2=18 parameters for each subject.

%% clear and path setup
clear all;close all;clc

rng('default'); % set the random seed to make the result reproducible.

postfix = 'RTmodel';

%% import one data
allData = load('taskDifficutlyOnIORExp1.mat');
subjList = unique(allData.Subject);
nSubj = length(subjList);

%% create data matrix for one subject
paramsfit = zeros(nSubj, 21);% column, 1-18,paramsters; 19,poslikelihood; 20, AIC; 21, BIC
paramsfit_all = cell(1,nSubj);

% use parallel computing toolbox, this command might depend on the Matlab
% version
parfor iSubj =1:nSubj  % loop subject
    iSubj
    % extract data and convert it into a matrix
    idx = find(allData.Subject==subjList(iSubj));
    data = {allData.ACC(idx) allData.RT(idx) allData.cf(idx) allData.cueing(idx) allData.difficulty(idx)};
    % data is a trial x five matrix. The column represents accuracy (1/0),
    % RT(ms), cf, cueing and difficulty
    % we replace string to number for three experimental variables
    data{3} = double(cell2mat(cellfun(@(x) strcmp(x,'on'), data{3}, 'UniformOutput',0))); % cf:1,on;0,off
    data{4} = double(cell2mat(cellfun(@(x) strcmp(x,'same'), data{4}, 'UniformOutput',0)));% cueing:1,same;0,diff
    data{5} = double(cell2mat(cellfun(@(x) strcmp(x,'hard'), data{5}, 'UniformOutput',0)));% difficulty:1,hard;0,easy
    data = cell2mat(data);
    
    % run the fitting
    out = fitRTmodel_optimize(data);
    [~,idx]=min(out(:,end-2)); % to get the minimum positive loglikelihood
    paramsfit_all{iSubj} = out;
    paramsfit(iSubj,:) = out(idx,:);
    %mybar(1:8,[paramsfit(3:10);paramsfit(11:18)])
end

t = fix(clock);
save(sprintf('%02d%02d%02d%02d%02d%02d_%s_%03dsubjs',t(1),t(2),t(3),t(4),t(5),t(6),postfix, nSubj));
