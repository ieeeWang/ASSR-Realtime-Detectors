clc; clear; close all;
% A demo using SPRT on 1 channel EEG

% add path for functions
utility_dir = '.\utility\';
addpath(genpath(utility_dir));

% load the FCz channel
load(['.\data\','FCz_100x.mat']) % [12s*100 trials]
EEGch = FCz;
fs = 2000;

%% compute F
F = []; % [100 trial * 999 Hz], not dB
for i=1:100  
    [Y_f1, f] = get_fftAmpSpec(EEGch(:,i), fs); 
    F(i,:) = get_spectrum_F(f, Y_f1, 12, 0.5);
end

%% choose ASSR freq to show
f_assr = [43, 80, 117];
F_tar = F(:,f_assr);

% check nan raw
N_nan = sum(isnan(F_tar(:,1)));
fprintf('there is %d nan raws in F matrix \n',N_nan);
if N_nan>0   
    % remove nan raws
    F_tar(isnan(F_tar(:,1)),:)=[];
    fprintf('removed... \n');
end

% check the Distribution of F
figure
nbins = 20;
for i=1:3
    histogram(F_tar(:,i), nbins);
    hold on
end
legend(num2str(f_assr'))
xlabel('F')
ylabel('Count')
title('Distribution of the power of generated f bins')

%% SPRT
n=1; % compute every n points
SNRdB_lower = -13;% dB
M = 12;
Utils = Utils_Bayesian;

tic
LLR = [];
SNR_ML_step = [];
for i=1:size(F_tar,2)
    [LLR(i,:), Nt, SNR_ML_step(i,:)] = Utils.SPRT_on_F_ML_SNR(F_tar(:,i), ...
        n, M, SNRdB_lower);
end
toc

%%
upper_lim = log(0.8/0.05);
lower_lim = -1.56;

% PLOT 3 ASSRs
legendtext = {'MF','2MF','3MF', 'upper bound','lower bound'};
marker1 = {'-o','-+','-x'};
legendtext2 = {'MF','2MF','3MF'};
fsize = 12;
LW=1;

figure
    subplot(211)    
    for i=1:3
        plot(Nt, LLR(i,:), 'LineWidth', LW); hold on
    end
    plot(Nt, upper_lim*ones(size(Nt)), '--.k');
    plot(Nt, lower_lim*ones(size(Nt)), '--k');
    grid on
    set(gca,'XTickLabel',[]) %removing numbers on axis without removing the grid
    ylim([-5 30])
    legend (legendtext, 'Location','northeast','NumColumns',1)
%     xlabel('Number of trials')
    ylabel('LLR', 'FontSize', fsize)
    title('SPRT with ML', 'FontSize', fsize)    
subplot(212) 
    plot(Nt, SNR_ML_step(1:3,:), 'LineWidth', LW);hold on
    grid on
%     legend (legendtext2, 'Location','southeast','NumColumns',1)
    ylim([-15 10])
    xlabel('# trials', 'FontSize', fsize)
    ylabel('SNR (dB)', 'FontSize', fsize)
    
    
%% get the step number where it surpass the threshold
thr = upper_lim;
Nassr = size(LLR,1)
step_NP = [];
for i=1:Nassr
    step_NP(i) = Utils.getStep_surpassThreshold(LLR(i,:), thr, 10); 
end
step_NP    
    
    
    
    
    
    
    
    
    
    
    