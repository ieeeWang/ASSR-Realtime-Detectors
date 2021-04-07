clc; clear; close all;
% 'Towards real-time detection of auditory steady-state responses: a
% comparative study' -- simulation of NP detector 

% add path for functions
utility_dir = '.\utility\';
addpath(genpath(utility_dir));

%% same with the paper M = 12;
SNRdB = linspace(-30,20,51); 
% from dB to power ratio ( SNR = signal power/noise power)
SNR = 10.^(SNRdB./10); 

Utils = Utils_Bayesian;
M = 12; % number of neighboring f bins
% compute simulation PD
muHat = []; muCI= []; PD= [];
muHat_log = []; muCI_log = [];
for i = 1:length(SNR)  
    [muHat(i), muCI(i), PD(i), crit_NP(i)]=Utils.NP_detector_SNR(SNR(i), M, 0);
    [muHat_log(i), muCI_log(i), ~, ~]=Utils.NP_detector_SNR(SNR(i), M, 1);
end

% compute theory PD
P_FA = 0.05;
PD_t = Utils.Pd_F_distribution(P_FA, M, SNR); 
    
%%
F_theory = SNR+1; % expected F score
fsize = 12;
figure
    subplot(211)
    errorbar(SNRdB, 10*log10(F_theory), muCI_log); hold on;
    plot(SNRdB, 10*log10(crit_NP),'-.');
    ylim([-5 23])
    grid on
    set(gca,'XTickLabel',[]) %removing numbers on axis without removing the grid lines
    legend('F (dB)','NP_c_r_i_t (\alpha = 0.05)', 'Location','northwest')
    ylabel('F (dB)', 'FontSize', fsize)
    title('NP detector', 'FontSize', fsize)
    
    subplot(212)
    plot(SNRdB, PD,'-o'); hold on;
    plot(SNRdB, PD_t);
    plot(SNRdB, 0.05*ones(size(SNRdB)),'--k');
    grid on
    legend({'simulation','theory', 'P_D = 0.05'}, 'Location','northwest')
    xlabel('SNR (dB)', 'FontSize', fsize)
    ylabel('P_D', 'FontSize', fsize)

%% compare multiple M values, i.e., number of neighboring f bins
M = [4, 6, 12, 24];
tic
muHat = []; muCI= []; PD_arr= [];
for m = 1:length(M)
    for i = 1:length(SNR)  
        [muHat(m,i), muCI(m,i), PD_arr(m,i), crit_NP(m,i)] = Utils.NP_detector_SNR(...
            SNR(i), M(m), 1);
    end
end
toc

% print the crit_NP under different M values
crit_NP1 = crit_NP(:,1)
crit_NP1_dB = 10*log10(crit_NP1)

% compute theory PD
P_FA = 0.05;
PD_t_array = [];
for i = 1:length(M)
    PD_t_array(i,:) = Utils.Pd_F_distribution(P_FA, M(i), SNR); 
end

%% 
% http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html#:~:text=Default%20Colors%20in%202D%20Graphs,-The%20default%20colors&text=In%20the%20past%2C%20each%20new,to%20manually%20adjust%20the%20colors.
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250];	
c4 = [0.4940, 0.1840, 0.5560];
mycolor = {c1,c2,c3,c4};% = default color in order;

legendtext1a = {'M_1=4','M_2=6','M_3=12','M_4=24'};
legendtext1b = {'NP_c_r_i_t_1 ','NP_c_r_i_t_2','NP_c_r_i_t_3','NP_c_r_i_t_4'};

legendtext2a = {'theory (M_1)','theory (M_2)','theory (M_3)','theory (M_4)'};
legendtext2b = {'simulation (M_1)','simulation (M_2)','simulation (M_3)','simulation (M_4)'};

FdB = 10*log10(F_theory);
FdB_arr = repelem(FdB,[4],[1]);

xrange = [-20, 20];
figure
    subplot(211)
    for i=1:4
        errorbar(SNRdB, FdB_arr(i,:), muCI(i,:)); hold on;
    end
    for i=1:4
        plot(SNRdB, 10*log10(crit_NP(i,:)),'-.', 'Color',mycolor{i}); hold on;
    end    
    ylim([-5 30])
    xlim(xrange)
    grid on
    set(gca,'XTickLabel',[]) %removing numbers on axis without removing the grid
    legend ([legendtext1a,legendtext1b],'NumColumns', 2,'Location','northwest')
    ylabel('F (dB)', 'FontSize', fsize)
    title('NP detector', 'FontSize', fsize)
    
    subplot(212)
    plot(SNRdB, PD_t_array); hold on
    plot(SNRdB, PD_arr(1,:), '*', 'Color',mycolor{1}); 
    plot(SNRdB, PD_arr(2,:), '.', 'Color',mycolor{2}); 
    plot(SNRdB, PD_arr(3,:), '+', 'Color',mycolor{3}); 
    plot(SNRdB, PD_arr(4,:), 'o', 'Color',mycolor{4}); 
    xlim(xrange)
    legend([legendtext2a, legendtext2b],'Location','northwest',...
        'NumColumns', 2)
    grid on
    xlabel('SNR (dB)', 'FontSize', fsize)
    ylabel('P_D', 'FontSize', fsize)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    