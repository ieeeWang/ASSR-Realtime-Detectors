clc; clear; close all;
% 'Towards real-time detection of auditory steady-state responses: a
% comparative study' -- simulation of SPRT detector 

% add path for functions
utility_dir = '.\utility\';
addpath(genpath(utility_dir));

%% a set of SNRs
M = 12;
SNRdB = [5, 0, -5, -10]; % dB
% SNRdB = [0, -5, -10, -20]; % dB
SNR_true = 10.^(SNRdB./10); % power ratio
N = 1000; %max No. of trials <<<<<<<<<<<<<<<<<
n = 10;% compute every 10 points 
SNRdB_lower = -13;% dB

tic
Utils = Utils_Bayesian;
n=10;% compute every 10 points
% implement SPRT on generated data from F distributions (F-test)
LR_1 = []; LR_0=[]; 
SNR_ML_step = [];
for i = 1:length(SNR_true)  
    % data under H1
    rng(10)
    Fscore_1 = ncfrnd(2, 2*M, 2*SNR_true(i), 1,N);
    [LR_1(i,:), Nt, SNR_ML_step(i,:)] = Utils.SPRT_on_F_ML_SNR(Fscore_1, n, M, SNRdB_lower);
end
% data under H0
rng(3)
Fscore_0 = frnd(2, 2*M, 1,N);
[LR_0, ~, SNR_ML_step(i+1,:)] = Utils.SPRT_on_F_ML_SNR(Fscore_0, n, M, SNRdB_lower);
toc

%% log plot 
legendtxt = {'5 dB', '0 dB', '-5 dB', '-10 dB', 'no signal', 'upper bound', 'lower bound'};
     
beta = 0.2; 
alpha = 0.05;
% alpha = 0.01;

upper_lim = log((1-beta)/alpha);
lower_lim = log(beta/(1-alpha));
fsize = 12;
LW=1;
figure
    subplot(211) 
    semilogx(Nt, LR_1, 'LineWidth', LW);hold on
    semilogx(Nt, LR_0,'LineWidth', LW);
    semilogx(Nt, upper_lim*ones(size(Nt)), '--.k');
    semilogx(Nt, lower_lim*ones(size(Nt)), '--k');
    grid on
    set(gca,'XTickLabel',[]) %removing numbers on axis without removing the grid
    ylim([-5 20])
    legend (legendtxt, 'Location','northwest','NumColumns',2)
    ylabel('LLR', 'FontSize', fsize)
    title('modified SPRT', 'FontSize', fsize) 
    
    subplot(212) 
    semilogx(Nt, SNR_ML_step, 'LineWidth', LW);hold on
    ylim([-15 10])
    grid on
    xlabel('# trials', 'FontSize', fsize)
    ylabel('SNR (dB)', 'FontSize', fsize)  
   
    
    
    
    
    
    
    