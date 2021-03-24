clc; clear; close all;
% NP vs SPRT in terms of Pd

%% load SPRT Pd
load(['.\data\','Pd_SPRT_SNR5_beta2.mat'])
SNRdB = [5, 0, -5, -10, -20]; % dB 
legendtxt = {'5 dB', '0 dB', '-5 dB', '-10 dB', '-20 dB'};


%% NP detector
% Ntrials = linspace(1,10001,10001); 
Ntrials = linspace(1,1001,1001); 

Utils = Utils_Bayesian;
tic
SNR_avg=[];Pd=[];
for i = 1:length(SNRdB)
    SNR_i = SNRdB(i);

    % compute SNR (dB) of the tempol averaging from N trials
    SNR_avg_tmp = 10*log10(Ntrials)+SNR_i;
    SNR_avg(:,i) = SNR_avg_tmp;
    % from dB to power ratio:
    SNR2 = 10.^(SNR_avg_tmp./10); % SNR2 = signal power/noise power

    P_FA = 0.05;
    M = 12;
    Pd(:,i) = Utils.Pd_F_distribution(P_FA, M, SNR2); 
end
toc
   

figure  
LW=1;
fsize = 12;
subplot(211)
    semilogx(Ntrials, Pd, 'LineWidth', LW)
    grid on
    set(gca,'XTickLabel',[]) %removing numbers on axis without removing the grid lines
    legend (legendtxt, 'Location','southeast','NumColumns',2)  
    xlim([10^1, 10^3])
    ylim([0 1.1])
    ylabel('P_D', 'FontSize', fsize)
    title('NP', 'FontSize', fsize)   
subplot(212)
    for i = 1:length(outputset)
        tmp = outputset{i};
        semilogx(tmp.Nt, tmp.Pd, 'LineWidth', LW); hold on
    end
    ylim([-0.1 1])
    grid on
    xlabel('# trials', 'FontSize', fsize)
    ylabel('P_D', 'FontSize', fsize)
    title('SPRT', 'FontSize', fsize)
    


    