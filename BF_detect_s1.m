clc; clear; close all;
% 'Towards real-time detection of auditory steady-state responses: a
% comparative study' -- simulation of BF detector 
% compute BF with thresholds (\eta=1, 3, 6)
% compare BF with NP

%% same with the paper M = 12;
SNRdB = linspace(-30,20,51); 
% from dB to power ratio:
SNR = 10.^(SNRdB./10); % SNR = signal power/noise power
F_theory = SNR+1; % expected F score
M = 12;

%% NP detector
Utils = Utils_Bayesian;

muHat = []; muCI= []; PD= [];
muHat_log = []; muCI_log = [];
for i = 1:length(SNR)  
    [muHat(i), muCI(i), PD(i), crit_NP(i)]=Utils.NP_detector_SNR(SNR(i), M, 0);
end

% compute theory PD
PD_t = Utils.Pd_F_distribution(0.05, M, SNR); 

%% Compute BF thresholds from PDFs of F  
[th_BF, p1] = Utils.getBayes_threshold_bf(SNR, M, 1);
[th_BF3, p3] = Utils.getBayes_threshold_bf(SNR, M, 3);
[th_BF6, p6] = Utils.getBayes_threshold_bf(SNR, M, 6);

figure
    subplot(211)
    plot(SNRdB, F_theory); hold on;
    plot(SNRdB, th_BF);
    plot(SNRdB, th_BF3);
    plot(SNRdB, th_BF6);
%     ylim([0 100])
    grid on
    legend('F','\theta_B_F(\eta=1)','\theta_B_F(\eta=3)', '\theta_B_F(\eta=6)',...
        'Location','northwest')
    xlabel('SNR (dB)')
    ylabel('F score')
    title('Bayes-factor detector')
    subplot(212)
    plot(SNRdB, 10*log10(F_theory)); hold on;
    plot(SNRdB, 10*log10(th_BF));
    plot(SNRdB, 10*log10(th_BF3));
    plot(SNRdB, 10*log10(th_BF6));
    grid on
    legend('F','\theta_B_F(\eta=1)','\theta_B_F(\eta=3)', '\theta_B_F(\eta=6)',...
        'Location','northwest')
    xlabel('SNR (dB)')
    ylabel('10*log_1_0(F)')
  
%% compute Pd and P_FA, PE of BF detector
output = Utils.perf_Bayes_detect(SNR, M, th_BF);
PD_b = output.Pd;
PE = output.PE;
P_FA = output.PFA;

output = Utils.perf_Bayes_detect(SNR, M, th_BF3);
PD_b3 = output.Pd;
PE3 = output.PE;
P_FA3 = output.PFA;

output = Utils.perf_Bayes_detect(SNR, M, th_BF6);
PD_b6 = output.Pd;
PE6 = output.PE;
P_FA6 = output.PFA;

%% compute confid (AC) from SNR
confid = Utils.getConfid(SNR, M);

%%
fsize = 12;
figure
    subplot(211)
    plot(SNRdB, 10*log10(F_theory)); hold on;
    plot(SNRdB, 10*log10(th_BF),'--');
    plot(SNRdB, 10*log10(th_BF3),'--.');  
    plot(SNRdB, 10*log10(th_BF6),'--*');  
    ylim([-3, 20])
    grid on
    legend('E(F)','\theta_F(\eta=1)','\theta_F(\eta=3)','\theta_F(\eta=6)',...
        'Location','northwest')
    ylabel('F (dB)', 'FontSize', fsize)
    set(gca,'XTickLabel',[]) %removing numbers on axis without removing the grid lines
    title('BF detector', 'FontSize', fsize)
    
    subplot(212)
    plot(SNRdB, confid); hold on;
    plot(SNRdB, 1-PE, '--'); 
    plot(SNRdB, 1-PE3, '--.');   
    plot(SNRdB, 1-PE6, '--*');
    ylim([0.5 1])
    grid on
    legend({'AC','PC (\eta=1)','PC (\eta=3)','PC (\eta=6)'},...
        'Location','northwest')
    xlabel('SNR (dB)', 'FontSize', fsize)
    ylabel('Probability', 'FontSize', fsize)


figure
    plot(SNRdB, confid); hold on;
    plot(SNRdB, PD_b, '--'); 
    plot(SNRdB, P_FA, '--'); 
    plot(SNRdB, PD_b3, '--.'); 
    plot(SNRdB, P_FA3, '--.'); 
    plot(SNRdB, PD_b6, '--*'); 
    plot(SNRdB, P_FA6, '--*'); 
    
    plot(SNRdB, PD_t,'-o'); 
    plot(SNRdB, 0.05*ones(size(SNRdB)),'--k');
    grid on
    legend({'AC (BF)','P_D (\eta=1)','P_F_A (\eta=1)',...
        'P_D (\eta=3)','P_F_A (\eta=3)',...
        'P_D (\eta=6)','P_F_A (\eta=6)',...
        'P_D (NP)', '\alpha=0.05 (NP)'},...
        'Location','northwest', 'NumColumns', 2)
    xlabel('SNR (dB)', 'FontSize', fsize)
    ylabel('Probability', 'FontSize', fsize)
    title('BF vs NP detectors', 'FontSize', fsize)






























    
    
    
    
    
    
    
    
    
    