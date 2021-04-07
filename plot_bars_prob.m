clc; clear; close all;
% plot fig for detecting all integer freqs

% add path for functions
utility_dir = '.\utility\';
addpath(genpath(utility_dir));

% load the ASSR frequencis
load(['.\data\','ASSR_4AMs.mat']) 
result_dir = '.\data\';

%% load outputset on 10 subj
load([result_dir,'output3_NP_BF_05.mat']) % alpha = 0.05
outputset1 = outputset;

ACmax=[];ACend=[]; 
PDmax =[]; PDend =[];
for j = 1:length(outputset1)
    PDmax(j,:) = outputset1{j}.PD_max;
    PDend(j,:) = outputset1{j}.PD_end;    
    ACmax(j,:) = outputset1{j}.AC_max;
    ACend(j,:) = outputset1{j}.AC_end;
end

% PDmax [9*100]
PD_O2 = PDmax(:,[1,2]); PD_O2 = reshape(PD_O2,1,[]);
PD_O4 = PDmax(:,[3:6]); PD_O4 = reshape(PD_O4,1,[]);
PD_O6 = PDmax(:,[7:12]); PD_O6 = reshape(PD_O6,1,[]);
PD_80 = PDmax(:,5);
PD_49 = PDmax(:,8);
PD_non = PDend(:,[13:end]); PD_non = reshape(PD_non,1,[]); % on 100th trials

PD_mean = [mean(PD_O2), mean(PD_O4), mean(PD_O6), mean(PD_80), mean(PD_49), mean(PD_non)];
PD_sd = [std(PD_O2,1), std(PD_O4,1), std(PD_O6,1), std(PD_80,1), std(PD_49,1), std(PD_non,1)];

% ACmax [9*100]
AC_O2 = ACmax(:,[1,2]); AC_O2 = reshape(AC_O2,1,[]);
AC_O4 = ACmax(:,[3:6]); AC_O4 = reshape(AC_O4,1,[]);
AC_O6 = ACmax(:,[7:12]); AC_O6 = reshape(AC_O6,1,[]);
AC_80 = ACmax(:,5);
AC_49 = ACmax(:,8);
AC_non = ACend(:,[13:end]); AC_non = reshape(AC_non,1,[]); % on 100th trials

AC_mean = [mean(AC_O2), mean(AC_O4), mean(AC_O6), mean(AC_80), mean(AC_49), mean(AC_non)];
AC_sd = [std(AC_O2,1), std(AC_O4,1), std(AC_O6,1), std(AC_80,1), std(AC_49,1), std(AC_non,1)];

% load PD (alpha = 0.01)
load([result_dir,'output3_NP_BF_01.mat']) % alpha = 0.01
outputset2 = outputset;
PDmax =[]; PDend =[];
for j = 1:length(outputset2)
    PDmax(j,:) = outputset2{j}.PD_max;
    PDend(j,:) = outputset2{j}.PD_end; 
end

% PDmax [9*100]
PD_O2 = PDmax(:,[1,2]); PD_O2 = reshape(PD_O2,1,[]);
PD_O4 = PDmax(:,[3:6]); PD_O4 = reshape(PD_O4,1,[]);
PD_O6 = PDmax(:,[7:12]); PD_O6 = reshape(PD_O6,1,[]);
PD_80 = PDmax(:,5);
PD_49 = PDmax(:,8);
PD_non = PDend(:,[13:end]); PD_non = reshape(PD_non,1,[]); % on 100th trials

PD_mean2 = [mean(PD_O2), mean(PD_O4), mean(PD_O6), mean(PD_80), mean(PD_49), mean(PD_non)];
PD_sd2 = [std(PD_O2,1), std(PD_O4,1), std(PD_O6,1), std(PD_80,1), std(PD_49,1), std(PD_non,1)];

%% plot PD, AC (m, std)
% tmp_avg=[PD_mean; AC_mean]; %[groups*bars]
% tmp_sd=[PD_sd; AC_sd];
tmp_avg=[PD_mean;PD_mean2; AC_mean]; %[groups*bars]
tmp_sd=[PD_sd;PD_sd2; AC_sd];

% multiple bars plot
figure
fsize = 12;
    numgroups = size(tmp_avg,1);  
    numbars = size(tmp_avg,2);  
    h = bar(tmp_avg);
    set(h,'BarWidth',0.8);    % The bars will now touch each other
    hold on;
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:numbars
          % Aligning error bar with individual bar
          x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  
          errorbar(x, tmp_avg(:,i), tmp_sd(:,i), 'k', 'linestyle', 'none');
    end
    ylim([0 1.2])
    legend({'Order 2','Order 4','Order 6','80 Hz','49 Hz','non-ASSR'},'Location','SouthEast')
%     set(gca,'Xtick',[1:2],'XTickLabel',{'P_D','AC'});
    set(gca,'Xtick',[1:3],'XTickLabel',{'P_D (\alpha = 0.05)','P_D (\alpha = 0.01)','AC'});
    set(gca,'ygrid','on')
%     xlabel('Groups','FontSize',fsize); 
    ylabel('Probability','FontSize',fsize);


