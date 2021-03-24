clc; clear; close all;
% A demo using BF vs. NP on 1 channel EEG

% add path for functions
utility_dir = '.\utility\';
addpath(genpath(utility_dir));

% load the FCz channel
load(['.\data\','FCz_100x.mat']) % [12s*100 trials]
EEGch = FCz;
fs = 2000;

%% compute F on accumulated trials
FdB = [];
for i=1:100
    tmp = EEGch(:,1:i);
    tmp_eeg = mean(tmp, 2, 'omitnan');
    
    [Y_f1, f] = get_fftAmpSpec(tmp_eeg, fs); 
    F = get_spectrum_F(f, Y_f1, 12, 0.5);
    FdB(i,:) = 10*log10(F);
end

%% choose ASSR freq
f_assr = [43, 80, 117];
FdB_tar = FdB(:,f_assr);
% move avarage along the columns
ws=5;
FdB_tar = movmean(FdB_tar,ws,1);

%% Part 1 (N-P):compute Pd
Utils = Utils_Bayesian;

% (1) FdB --> SNR estimate
SNR_tar = Utils.SNRdB_Nt2SNR_1Tr(FdB_tar, 1);

% (2) compute Pd from SNR
P_FA = 0.05;  M = 12;
P_D = [];
for i=1:size(SNR_tar,2)
    P_D(:,i) = Utils.Pd_F_distribution(P_FA, M, SNR_tar(:,i)); 
end

%% plot
% computed NP critical values (M=12), refer to \NP_detector_s1.m   
sig_SNR_05 = 5.318;% dB
sig_SNR_01 = 7.492;% dB
legendtext = {'MF','2MF','3MF'};
marker1 = {'-o','-+','-x','-.'};
% marker2 = {'--*','--+','--x','--.'};
fsize = 12;
LW=1;

figure
subplot(211)
    Tmax = 100;
    for i=1:3
        plot(FdB_tar(:,i),'LineWidth', LW); hold on
    end
    plot((1:Tmax),sig_SNR_05*ones(Tmax,1),'k--','LineWidth', 1)
    grid on
    set(gca,'XTickLabel',[]) %removing numbers on axis without removing the grid
    ylim([-20 30])
    ylabel('F (dB)', 'FontSize', fsize) 
    title('NP', 'FontSize', fsize)   
    legend([legendtext, 'NP_c_r_i_t'], 'Location','southeast')

subplot(212)
    for i=1:3
        plot(P_D(:,i),'LineWidth', LW); hold on
    end
    grid on
    ylim([0 1.1])
    xlabel('# trials', 'FontSize', fsize);
    ylabel('P_D', 'FontSize', fsize)    
%     legend(legendtext, 'Location','southeast')
    
    
%% Part 2: Bayes-factor
tic
% (1) compute B-F thresholds from SNR
M = 12;
th_BF = [];
for i=1:size(SNR_tar,2)
    [th_BF(:,i), ~] = Utils.getBayes_threshold(SNR_tar(:,i), M);% bf=1
    % Same with above
%     [th_BF(:,i), ~] = Utils.getBayes_threshold_bf(SNR_tar(:,i), M, 1);
    [th_BF3(:,i), ~] = Utils.getBayes_threshold_bf(SNR_tar(:,i), M, 3); % bf=3
end
toc
% (2) compute AC from SNR (estimate)
AC = [];
for i=1:size(SNR_tar,2)
    AC(:,i) = Utils.getConfid(SNR_tar(:,i), M);
end

%% (3) compute PC from SNR (estimate)
PC = [];
for i=1:size(SNR_tar,2)
    output = Utils.perf_Bayes_detect(SNR_tar(:,i), M, th_BF(:,i));
    PC(:,i) = 1- output.PE;
end

PC3 = [];
for i=1:size(SNR_tar,2)
    output = Utils.perf_Bayes_detect(SNR_tar(:,i), M, th_BF3(:,i));
    PC3(:,i) = 1- output.PE;
end

%% plot
% http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html#:~:text=Default%20Colors%20in%202D%20Graphs,-The%20default%20colors&text=In%20the%20past%2C%20each%20new,to%20manually%20adjust%20the%20colors.
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250];	
c4 = [0.4940, 0.1840, 0.5560];
mycolor = {c1,c2,c3,c4};% = default color in order;

legendtext1 = {'F (MF)','F (2MF)','F (3MF)', '\theta_F (MF)','\theta_F (2MF)','\theta_F (3MF)'};
legendtext2 = {'AC (MF)','AC (2MF)','AC (3MF)', 'PC (MF)','PC (2MF)','PC (3MF)'};

figure
    subplot(211)
    Tmax = 100;
    for i=1:3
        %plot(FdB_tar(:,i), marker1{i},'Color',mycolor{i}); hold on
        plot(FdB_tar(:,i),'Color',mycolor{i},'LineWidth', LW); hold on
    end
    for i=1:3
        plot(10*log10(th_BF(:,i)), '--', 'Color',mycolor{i});
    end
    grid on
    set(gca,'XTickLabel',[]) %removing numbers on axis without removing the grid
    legend (legendtext1,'NumColumns', 2, 'Location','southeast')
    ylim([-20 30])
    ylabel('F (dB)', 'FontSize', fsize) 
    title(['BF (\eta=1)'], 'FontSize', fsize)   
    
    subplot(212)    
    for i=1:3
        plot(AC(:,i), 'Color',mycolor{i},'LineWidth', LW); hold on
    end
    for i=1:3
        plot(PC(:,i),'--', 'Color',mycolor{i}); hold on
    end    
    plot((1:Tmax),0.5*ones(Tmax,1),'k--','LineWidth', 1)
    grid on
    ylim([0.4 1.1])
    legend(legendtext2, 'Location','southeast', 'NumColumns', 2)
    xlabel('# trials', 'FontSize', fsize);
    ylabel('Probability', 'FontSize', fsize)
    
%%    
figure
    subplot(211)
    Tmax = 100;
    for i=1:3
        %plot(FdB_tar(:,i), marker1{i},'Color',mycolor{i}); hold on
        plot(FdB_tar(:,i),'Color',mycolor{i},'LineWidth', LW); hold on
    end
    for i=1:3
        plot(10*log10(th_BF3(:,i)), '--', 'Color',mycolor{i});
    end
    grid on
    set(gca,'XTickLabel',[]) %removing numbers on axis without removing the grid
    legend (legendtext1,'NumColumns',2, 'Location','southeast')
    ylim([-20 30])
    ylabel('F (dB)', 'FontSize', fsize) 
    title(['BF (\eta=3)'], 'FontSize', fsize)   
    
    subplot(212)    
    for i=1:3
        plot(AC(:,i), 'Color',mycolor{i},'LineWidth', LW); hold on
    end
    for i=1:3
        plot(PC3(:,i),'--', 'Color',mycolor{i}); hold on
    end    
    plot((1:Tmax),0.75*ones(Tmax,1),'k--','LineWidth', 1)
    grid on
    ylim([0.4 1.1])
    legend(legendtext2, 'Location','southeast', 'NumColumns', 2)
    xlabel('# trials', 'FontSize', fsize);
    ylabel('Probability', 'FontSize', fsize)

%% print max F and corresponding PD AC
[FdB_tar_max, I]= max(FdB_tar,[],1);
FdB_tar_max
I

P_D_max = [];
for i=1:length(I)
    P_D_max(i) = P_D(I(i),i);
end
P_D_max

AC_max = [];
for i=1:length(I)
    AC_max(i) = AC(I(i),i);
end
AC_max




















    