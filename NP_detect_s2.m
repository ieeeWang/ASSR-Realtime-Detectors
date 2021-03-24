clc; clear; close all;
% 'Towards real-time detection of auditory steady-state responses: a
% comparative study' -- simulation of NP detector 
% plot contour of Pd (i.e., function of trial number and SNR)
%%
% Ntrials = linspace(1,10001,10001); % 386 s
Ntrials = linspace(1,1001,1001); % 34 s

% SNRdB = linspace(-30,10,5); % every 10 dB
% SNRdB = linspace(-30,10,9); % every 5 dB
SNRdB = linspace(-30,10,41); % every 1 dB

tic
SNR_avg=[];Pd=[];
for i = 1:length(SNRdB)
    SNR_i = SNRdB(i);

    % compute SNR (dB) of the tempol averaging from N trials
    SNR_avg_tmp = 10*log10(Ntrials)+SNR_i;
    SNR_avg(:,i) = SNR_avg_tmp;
    % from dB to power ratio:
    SNR2 = 10.^(SNR_avg_tmp./10);

    Utils = Utils_Bayesian;
    P_FA = 0.05;  M = 12;
    Pd(:,i) = Utils.Pd_F_distribution(P_FA, M, SNR2); 
end
toc

%%
fsize = 12;
dbtxt = {'-30 dB','-20 dB','-10 dB','0 dB','10 dB'}';

figure
subplot(211)
%     plot(Ntrials, SNR_avg)
    semilogx(Ntrials, SNR_avg)
    grid on
    legend (dbtxt, 'Location','southeast')    
    ylabel('SNR (dB)', 'FontSize', fsize)
    set(gca,'XTickLabel',[]) %removing numbers on axis without removing the grid lines
%     title('SNR as a function of # trials')
subplot(212)
%     plot(Ntrials, Pd)
    semilogx(Ntrials, Pd)
    grid on
    %legend ([cellstr(num2str(SNRdB'))], 'Location','southeast')  
    xlabel('# trials', 'FontSize', fsize)
    ylabel('P_D', 'FontSize', fsize)
%     title('P_D as a function of # trials')

%% 3-D plot & contour
fsize = 13;
h=figure;
    A = axes;
    % contour(Ntrials,SNRdB,Pd', 'ShowText','on');
    contourf(Ntrials,SNRdB,Pd','ShowText','on');
    colormap default    % change color map
    set(A, 'xScale', 'log')
    xlabel('# trials', 'FontSize', fsize)
    ylabel('SNR (dB)', 'FontSize', fsize)
    title('Contour of P_D', 'FontSize', fsize) 

%% compute Minimum single-trial SNRs (dB) needed for achieving a Pd
% make sure: SNRdB = linspace(-30,10,41); % every 1 dB
Tmax = 100;
Pd_Tmax = Pd(Tmax,:);

loc = [];
Pd_target = [0.5, 0.6, 0.7, 0.8, 0.9, 1];
for i = 1:length(Pd_target)
    tmp = Pd_Tmax - Pd_target(i);
    [M,I] = min(abs(tmp));
    loc(i)= I;
end
SNRmin = SNRdB(loc)

