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
load([result_dir,'output2_NP_BF_01.mat'])
outputset1 = outputset;

FdBmax=[];FdBend=[];
ACmax=[];ACend=[];
for j = 1:length(outputset1)
    FdBmax(j,:) = outputset1{j}.FdB_max;
    FdBend(j,:) = outputset1{j}.FdB_end;   
    ACmax(j,:) = outputset1{j}.AC_max;
    ACend(j,:) = outputset1{j}.AC_end;
end

load([result_dir,'output2_SPRT_01.mat'])
outputset2 = outputset;
LLRmax=[];LLRend=[];
for j = 1:length(outputset2)
    LLRmax(j,:) = outputset2{j}.LLR_max;
    LLRend(j,:) = outputset2{j}.LLR_end;
end

load([result_dir,'output2_T2.mat'])
outputset3 = outputset;
T2=[];
for j = 1:length(outputset3)
    T2(j,:) = outputset3{j}.T2;
end

    
%% Default Colors in 2D Graphs (default colors in MATLAB changed in R2014b version)
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250];	
c4 = [0.4940, 0.1840, 0.5560];
mycolor = {c1,c2,c3,c4};% default color order;


%% plot F & AC 
% NP crit
% computed critical values, refer to \NP_detect_s1.m   
sig_SNR_05 = 5.32;% dB
sig_SNR_01 = 7.49;% dB
% AC thresholds
BF3 = 0.75; % BF=3
BF6 = 6/7; % BF=6
% LLR thresholds
beta = 0.2; % PD = 1-beta
alpha = 0.01;
upper_lim = log((1-beta)/alpha);
lower_lim = log(beta/(1-alpha));
% significance level of T^2, refer to reference paper
sig_T2_05 = 6.1; % p = 0.05
sig_T2_01 = 9.6; % p = 0.01


figure  
    sn = 12;
    fmax = 150;
subplot(411)
    stem(f_target1, mean(FdBmax(:,1:12),1),'*', 'Color',mycolor{1}); hold on % ASSRS
    stem(f_target0, mean(FdBend(:,13:end),1), 'Color',mycolor{1}) % non-ASSRS 
    plot((1:fmax),sig_SNR_01*ones(fmax,1),'r--','LineWidth', 1)
    plot((1:fmax),sig_SNR_05*ones(fmax,1),'k--','LineWidth', 1) 
    grid on 
    set(gca,'XTickLabel',[]) %removing numbers on axis without removing the grid lines
    ylim([0 25])
    legend ('ASSRs','non-ASSRs', 'NP_c_r_i_t (\alpha=0.01)', ...
        'NP_c_r_i_t (\alpha=0.05)')
    ylabel('F (dB)','Fontsize',sn)
    title('NP','Fontsize',sn)
subplot(412)
    plot((1:fmax),BF6*ones(fmax,1),'r--','LineWidth', 1); hold on
    plot((1:fmax),BF3*ones(fmax,1),'k--','LineWidth', 1)
    stem(f_target1, mean(ACmax(:,1:12),1),'*','Color',mycolor{2}) % ASSRS
    stem(f_target0, mean(ACend(:,13:end),1),'Color',mycolor{2}) % non-ASSRS        
    grid on 
    set(gca,'XTickLabel',[]) 
    ylim([0.5 1.2])
    legend ('\theta_A_C (\eta=6)', '\theta_A_C (\eta=3)')
    ylabel('AC','Fontsize',sn)
    title('BF','Fontsize',sn)
subplot(413)
    plot((1:fmax),upper_lim*ones(fmax,1),'r--','LineWidth', 1); hold on
    plot((1:fmax),lower_lim*ones(fmax,1),'k--','LineWidth', 1)
    stem(f_target1, mean(LLRmax(:,1:12),1),'*','Color',mycolor{3}) % ASSRS
    stem(f_target0, mean(LLRend(:,13:end),1),'Color',mycolor{3}) % non-ASSRS        
    grid on 
    set(gca,'XTickLabel',[])
    ylim([-5 15])
    legend ('upper bound', 'lower bound')
    ylabel('LLR','Fontsize',sn)
    title('modified SPRT','Fontsize',sn)    
subplot(414)  
    plot((1:fmax),sig_T2_01*ones(fmax,1),'r--','LineWidth', 1); hold on
    plot((1:fmax),sig_T2_05*ones(fmax,1),'k--','LineWidth', 1)
    stem(f_target1, mean(T2(:,1:12),1),'*','Color',mycolor{4}) % ASSRS
    stem(f_target0, mean(T2(:,13:end),1),'Color',mycolor{4}) % non-ASSRS   
    ylim([0 20]) 
    grid on 
    legend ('sig. level p=0.01','sig. level p=0.05')
    ylabel('T^2','Fontsize',sn)   
    xlabel('Frequency (Hz)','Fontsize',sn)
    title('Hotteling T2','Fontsize',sn)
    xAxis_number =  [0, 37, 43, 50, 80, 100, 150];
    xAxis_text = {0, 'f_1', 'f_2', 50, 80, 100, 150};
    xticks(xAxis_number);xticklabels(xAxis_text)
    
    
%% plot P-R curves: on accumulated values across subjs
label_bi = zeros([1, 100]); label_bi([1:12]) = 1;

% remove two line noise freq @ 50 and 100 Hz, does not affect anything!
f1 = find(f_target == 50);
f2 = find(f_target == 100);
% LLRend(:,[f1,f2])=0;

N = 9; % subj
x1=[]; y1=[];
for i=1:N
    tmp = FdBend(i,:);
    tmp([1:12])=FdBmax(i,[1:12]); % ASSRs
    x1 = [x1, tmp];
    y1 = [y1, label_bi];
end

x2=[]; 
for i=1:N
    tmp = ACend(i,:);
    tmp([1:12])=ACmax(i,[1:12]); % ASSRs
    x2 = [x2, tmp];
end

x3=[];
for i=1:N
    tmp = LLRend(i,:);
    tmp([1:12])=LLRmax(i,[1:12]); % ASSRs
    x3 = [x3, tmp];
end

x4=[]; % T2
for i=1:N
    x4 = [x4, T2(i,:)];
end


plt = 0;
ylim_range = [0 1];
r1 = ROC2plot(x1, y1, plt,0,ylim_range);
r2 = ROC2plot(x2, y1, plt,0,ylim_range);
r3 = ROC2plot(x3, y1, plt,0,ylim_range);
r4 = ROC2plot(x4, y1, plt,0,ylim_range);

figure
    fn = 14;
% plot multiple PR curves in one figure
subplot(121)
    plot(r1.FPrate, r1.sensitivity, 'Color',mycolor{1},'LineWidth',2); hold on;
    plot(r2.FPrate, r2.sensitivity, 'Color',mycolor{2},'LineWidth',2);
    plot(r3.FPrate, r3.sensitivity, 'Color',mycolor{3},'LineWidth',2);
    plot(r4.FPrate, r4.sensitivity, 'Color',mycolor{4},'LineWidth',2);
    grid on
    xlim([0 1])
    ylim(ylim_range)
    ylabel('Sensitivity','FontSize',fn);
    xlabel('1- specificity','FontSize',fn);
    legend(['NP: AUC = ',sprintf('%.2f', r1.auc)],...
        ['BF: AUC = ',sprintf('%.2f', r2.auc)],...
        ['SPRT: AUC = ',sprintf('%.2f', r3.auc)],...
        ['T^2: AUC = ',sprintf('%.2f', r4.auc)],...
        'Location','southwest')
    title(['ROC Curves'],'FontSize',fn);

    
ylim_range = [0 1];
r1 = PR2plot(x1, y1, plt,0,ylim_range);
r2 = PR2plot(x2, y1, plt,0,ylim_range);
r3 = PR2plot(x3, y1, plt,0,ylim_range);
r4 = PR2plot(x4, y1, plt,0,ylim_range);

% plot multiple PR curves in one figure
subplot(122)
    plot(r1.recall, r1.precision,'Color',mycolor{1},'LineWidth',2); hold on;
    plot(r2.recall, r2.precision,'Color',mycolor{2},'LineWidth',2);
    plot(r3.recall, r3.precision,'Color',mycolor{3},'LineWidth',2);
    plot(r4.recall, r4.precision,'Color',mycolor{4},'LineWidth',2);
    grid on
    xlim([0 1])
    ylim(ylim_range)
    xlabel('Recall (sensitivity)','FontSize',fn);
    ylabel('Precision (PPV)','FontSize',fn);
    legend([sprintf('%.2f', r1.auc)],...
        [sprintf('%.2f', r2.auc)],...
        [sprintf('%.2f', r3.auc)],...
        [sprintf('%.2f', r4.auc)],...
        'Location','southwest')
    title(['P-R Curves'],'FontSize',fn);
  

    
    
