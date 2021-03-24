clc; clear; close all;
% plot fig for detecting all integer freqs

% add path for functions
utility_dir = '.\utility\';
addpath(genpath(utility_dir));

% load the ASSR frequencis
load(['.\data\','ASSR_4AMs.mat']) 
result_dir = '.\data\';

%% load ASSR f
f_target1 = [37, 43, 6, 74, 80, 86, 31, 49, 111, 117, 123, 129]; % Fc = 500 Hz
f_MB = [ASSR_4AMs.MB2; ASSR_4AMs.MB3];
f_MB = unique(f_MB);
f_BB = [ASSR_4AMs.BB2; ASSR_4AMs.BB3];

f_fullset = [2:149];
tmpset = setdiff(f_fullset, f_MB); % setA - setB
f_target0 = setdiff(tmpset, f_BB);
f_target = [f_target1, f_target0]; % integer freqs: (non-ASSR + ASSR)

%% load outputset on 10 subj
% load([result_dir,'output3_NP_BF_01.mat']) % alpha = 0.01
load([result_dir,'output3_NP_BF_05.mat']) % alpha = 0.05
outputset1 = outputset;

FdBmax=[];FdBend=[];
ACmax=[];ACend=[];
for j = 1:length(outputset1)
    FdBmax(j,:) = outputset1{j}.FdB_max;
    FdBend(j,:) = outputset1{j}.FdB_end;   
    ACmax(j,:) = outputset1{j}.AC_max;
    ACend(j,:) = outputset1{j}.AC_end;
end

% load([result_dir,'output3_SPRT_01.mat']) % alpha = 0.01
load([result_dir,'output3_SPRT_05.mat']) % alpha = 0.05
outputset2 = outputset;
LLRmax=[];LLRend=[];
for j = 1:length(outputset2)
    LLRmax(j,:) = outputset2{j}.LLR_max;
    LLRend(j,:) = outputset2{j}.LLR_end;
end

%% print for the latex table
for j = 1:length(outputset1)
    tmp = outputset1{j};
    fprintf('Subj#%d =====================================\n', j) 
    step_maxF = tmp.step_maxF;
    FdB = tmp.FdB_max;
    N =12; % ASSRs
    for i= 1: N 
        fprintf('&%3.0f(%d) ', FdB(i), step_maxF(i)) 
    end
    fprintf('FdB \n')    
    
    step_NP = tmp.step_NP;
    PD = tmp.PD_max;
    for i= 1: N
        fprintf('&%3.2f(%d) ', PD(i), step_NP(i)) 
    end
    fprintf('NP \n')

    step_BF = tmp.step_BF;
    AC = tmp.AC_max;
    for i= 1:N
        fprintf('&%3.2f(%d) ', AC(i), step_BF(i)) 
    end
    fprintf('BF \n')
    
    tmp = outputset2{j};    
    for i= 1: N
        fprintf('& $Z_1$(%d) ', tmp.step_SPRT2(i)) 
    end
    fprintf('SPRT \n')
end


%% print perf. table
Utils = Utils_Bayesian;
% NP
for j = 1:length(outputset1) 
    tmp = outputset1{j}; 
    step_vec = tmp.step_NP2;
    % Get perf.
    confMatrix = Utils.stepVec2confMatrix(step_vec, length(f_target1));
    perf = Utils.getPefconfMatrix(confMatrix);
    recall(j) = perf.recall;
    precision(j) = perf.precision;
    specif(j) = perf.specif;
end
fprintf('==============NP: ==================\n') 
fprintf('recall mean: %3.2f, std: %3.2f \n', mean(recall), std(recall))
fprintf('precision mean: %3.2f, std: %3.2f \n', mean(precision), std(precision))
fprintf('specif mean: %3.2f, std: %3.2f \n', mean(specif), std(specif))

% BF
for j = 1:length(outputset1) 
    tmp = outputset1{j}; 
    step_vec = tmp.step_BF2;
    % Get perf.
    confMatrix = Utils.stepVec2confMatrix(step_vec, length(f_target1));
    perf = Utils.getPefconfMatrix(confMatrix);
    recall(j) = perf.recall;
    precision(j) = perf.precision;
    specif(j) = perf.specif;
end
fprintf('==============BF: ==================\n') 
fprintf('recall mean: %3.2f, std: %3.2f \n', mean(recall), std(recall))
fprintf('precision mean: %3.2f, std: %3.2f \n', mean(precision), std(precision))
fprintf('specif mean: %3.2f, std: %3.2f \n', mean(specif), std(specif))

% SPRT
for j = 1:length(outputset2)
    tmp = outputset2{j}; 
    step_vec = tmp.step_SPRT2;
    % Get perf.
    confMatrix = Utils.stepVec2confMatrix(step_vec, length(f_target1));
    perf = Utils.getPefconfMatrix(confMatrix);
    recall(j) = perf.recall;
    precision(j) = perf.precision;
    specif(j) = perf.specif;
end
fprintf('============== SPRT: ==================\n') 
fprintf('recall mean: %3.2f, std: %3.2f \n', mean(recall), std(recall))
fprintf('precision mean: %3.2f, std: %3.2f \n', mean(precision), std(precision))
fprintf('specif mean: %3.2f, std: %3.2f \n', mean(specif), std(specif))


%% for box plot 
F12_NP = []; F12_BF = []; F12_SPRT = []; 
N =12; % ASSRs
for j = 1:length(outputset1)
    tmp = outputset1{j};
    F12_NP(j,:) = tmp.step_NP(1:N);
    F12_BF(j,:) = tmp.step_BF(1:N);
    F12_SPRT(j,:) = outputset2{j}.step_SPRT2(1:N);
end

n_sig = [];
n_sig(1,:) = sum(~isnan(F12_NP));
n_sig(2,:) = sum(~isnan(F12_BF));
n_sig(3,:) = sum(~isnan(F12_SPRT));
n_sig(4,:) = f_target1;

%% box plot 1
MyUtilis = Utils_Bayesian; %
steps_O3_NP = MyUtilis.getsteps_orders(F12_NP); 
steps_O3_BF = MyUtilis.getsteps_orders(F12_BF); 
steps_O3_SPRT = MyUtilis.getsteps_orders(F12_SPRT); 

Y = cat(3, steps_O3_NP(:,[1:3]), steps_O3_BF(:,[1:3]), steps_O3_SPRT(:,[1:3]));
% adjust dim to: samples*Xgroups*LABELS, e.g.,Y = permute(Y, [1,3,2]);

X = [1 2 3];
figure
    fsize = 12;
    iosr.statistics.boxPlot(X,Y,...
        'symbolColor','k',...
        'medianColor','k',...
        'symbolMarker',{'o','o','o'},...
        'boxcolor',{[.25 .25 .25],[.5 .5 .5],[.75 .75 .75]},... %'boxcolor','auto',...
         'groupLabels',{'NP','BF','SPRT'},... % LABELS
         'showLegend',true,...
          'showScatter',true);
    set(gca,'XTickLabel',{'Order 2','Order 4','Order 6'}, 'FontSize', fsize)
    ylabel('# trials', 'FontSize', fsize) 
    title('Trials needed for ASSR detection', 'FontSize', fsize) 
    box on
    grid on

%% add 80 Hz (index = 5) and 49 Hz (i = 8)
Y = cat(3, steps_O3_NP, steps_O3_BF, steps_O3_SPRT);
X = [1 2 3 4 5];
figure
    fsize = 12;
    iosr.statistics.boxPlot(X,Y,...
        'symbolColor','k',...
        'medianColor','k',...
        'symbolMarker',{'o','o','o'},...
        'boxcolor',{[.25 .25 .25],[.5 .5 .5],[.75 .75 .75]},... %'boxcolor','auto',...
         'groupLabels',{'NP','BF','SPRT'},... % LABELS
         'showLegend',true,...
          'showScatter',true);
    set(gca,'XTickLabel',{'Order 2','Order 4','Order 6','80 Hz','49 Hz'}, 'FontSize', fsize)
    ylabel('# trials', 'FontSize', fsize) 
    title('Accumulated trials needed for detecting ASSRs', 'FontSize', fsize) 
    box on
    grid on











