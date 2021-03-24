clc; clear; close all;
% compute Pd of SPRT and save it for comparison with the NP detector

%% (1) generate data
N_MC = 1e2; % ~67 s 
% N_MC = 1e4; % ~20 mins

% max sample size of data we observe
Ns_max = 1000;% max trials
% the true SNR, to be estimate by ML estimator
SNRdB_true = [5, 0, -5, -10, -20]; % dB 

%%
tic
Utils = Utils_Bayesian;
n = 10;% compute every 10 points
M = 12; % n of neighbering fins

outputset={};
for i =1:length(SNRdB_true)
    outputset{i} = Utils.get_Pd_SPRT(SNRdB_true(i), Ns_max, N_MC, M, 1);
    toc
end
toc

%% save for plot
% save('Pd_SPRT_SNR5_beta2.mat', 'outputset')

%%
upper_lim = 2.77;
lower_lim = -1.56;
    
% upper_lim1 = log(0.5/0.05);
% upper_lim2 = log(0.8/0.05);

figure
    subplot(211)
    for i = 1:length(outputset)
        tmp = outputset{i};
        plot(tmp.Nt, tmp.Pd); hold on
    end
    legend(cellstr(num2str(SNRdB_true')))
    grid on
    ylim([-0.1 1])
    xlabel('# trials')
    ylabel('P_D')
    subplot(212)
    for i = 1:length(outputset)
        tmp = outputset{i};
        semilogx(tmp.Nt, tmp.Pd); hold on
    end
    legend(cellstr(num2str(SNRdB_true')))
    ylim([-0.1 1])
    grid on
    xlabel('# trials')
    ylabel('P_D')
    