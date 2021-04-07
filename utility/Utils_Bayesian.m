% helper functions for implementing Bayesian
% Lei Wang, ieeewangl@gmail.com
% Oct. 2020 - Jan. 2021 @Radboud Uni.

function Utils_Bayesian = Utils_Bayesian

    Utils_Bayesian.NP_detector_SNR = @NP_detector_SNR;
    Utils_Bayesian.Pd_F_distribution = @Pd_F_distribution;
    Utils_Bayesian.MAP_estimate_SNRofonetrial= @MAP_estimate_SNRofonetrial;
    Utils_Bayesian.SPRT_on_Fscore= @SPRT_on_Fscore;
    Utils_Bayesian.SPRT_on_F_ML_SNR= @SPRT_on_F_ML_SNR;
    Utils_Bayesian.get_Pd_SPRT= @get_Pd_SPRT;
    Utils_Bayesian.SNRdB_Nt2SNR_1Tr= @SNRdB_Nt2SNR_1Tr;
    Utils_Bayesian.find_intersect= @find_intersect;
    Utils_Bayesian.getBayes_threshold= @getBayes_threshold;
    Utils_Bayesian.getBayes_threshold_bf= @getBayes_threshold_bf;
    
    Utils_Bayesian.getConfid= @getConfid;
    Utils_Bayesian.perf_Bayes_detect= @perf_Bayes_detect;
    
    Utils_Bayesian.getStep_surpassThreshold= @getStep_surpassThreshold;
    Utils_Bayesian.step100th_surpassThreshold= @step100th_surpassThreshold;
    
    Utils_Bayesian.stepVec2confMatrix= @stepVec2confMatrix;
    Utils_Bayesian.pp2confMatrix = @pp2confMatrix;
    Utils_Bayesian.getPefconfMatrix= @getPefconfMatrix;
    Utils_Bayesian.getfft_compx= @getfft_compx;
    
    Utils_Bayesian.getsteps_orders= @getsteps_orders;
end



function output = getsteps_orders(tmp) 
% get the steps for 2, 4 and 6 orders
MF1 = [tmp(:,1); tmp(:,2)];
MF2 = [tmp(:,3); tmp(:,4); tmp(:,5); tmp(:,6)];
MF3 = [];
for i = 7:12
    MF3 = [MF3; tmp(:,i)];
end
L_max = length(MF3);
output = [padding_nan(MF1, L_max), padding_nan(MF2, L_max), MF3];

% add 80 Hz (index = 5) and 49 Hz (i = 8)
output = [output, padding_nan(tmp(:,5), L_max), padding_nan(tmp(:,8),L_max)];
end


function C1 = padding_nan(loc, L_max)
C1 = loc(:);
C1 = [C1; nan(L_max-length(loc),1)];
end


function [P_compx] = getfft_compx(X, fs, N) 
% get the complex DFT at an integer frequency N, (0<N<fs/2)

% f resolusion =  1/t, t is duration (sec) of X
% DFT is computed without setting a NFFT, see explanation in My_project\test\FFT_with_out_NFFT.m
% OUTPUT
% P_compx - complex value at the integer frequency N (Hz)


Y = fft(X); % same with Y = fft(X, lenth(X)); to ensure the maximum f resolusion
L = length(X); % f resolusion = fs/L = 1/t
Y_comp = Y/L;
P1 = Y_comp(1 : (L/2+1)); % first half including the median point
P1(2:end-1) = 2*P1(2:end-1);
f = fs/2*(0:(2/L):1); % same with f = Fs*(0:(L/2))/L;

% P_compx = P1((f==N)); % it causes error when 99 != 99.0000 
P_compx = P1(12*N+1); % use this when f resolution = 1/12 Hz
end


function confMatrix = pp2confMatrix(stepVec, firstN, thr)
% stepVec: firstN items are positive claass, negtive for remaining
confMatrix = {};
% the real positive vec:
TP_vec = stepVec(1:firstN);
TP =  sum(TP_vec < thr);
FN = length(TP_vec) - TP;
% the real negtive vec:
TN_vec = stepVec(firstN+1:end);
TN = sum(TN_vec>thr);
FP = length(TN_vec) - TN;

confMatrix.TP = TP;
confMatrix.FN = FN;
confMatrix.TN = TN;
confMatrix.FP = FP;

end


function confMatrix = stepVec2confMatrix(stepVec, firstN)
% stepVec: firstN items are positive claass, negtive for remaining
confMatrix = {};
% the real positive vec:
TP_vec = stepVec(1:firstN);
TP =  firstN - sum(isnan(TP_vec));
FN = sum(isnan(TP_vec));
% the real negtive vec:
TN_vec = stepVec(firstN+1:end);
TN = sum(isnan(TN_vec));
FP = length(TN_vec) - TN;

confMatrix.TP = TP;
confMatrix.FN = FN;
confMatrix.TN = TN;
confMatrix.FP = FP;

end


function perf = getPefconfMatrix(confMatrix)

TP= confMatrix.TP;
FN= confMatrix.FN;
TN= confMatrix.TN;
FP= confMatrix.FP;

perf = {};
perf.recall = TP/(TP+FN);
perf.precision = TP/(TP+FP);
perf.specif = TN/(TN+FP);

end


function step_sig = step100th_surpassThreshold(vec, thr)
% check if 100th F(dB) surpass threshold

if vec(end)>thr 
    step_sig = 100;
else
    step_sig = nan; % no surpass
end
end


function step_sig = getStep_surpassThreshold(vec, thr, minim_step)
% get the step number where it surpass the threshold, thr
% minim_step = 0 by default

step_sig = nan; % no surpass
for i = 1:length(vec)
    if vec(i)>thr && (i>minim_step)
        step_sig = i;
        break
    end
end

if ~isnan(step_sig)
    if step_sig<minim_step
        step_sig = minim_step;
    end
end
end


function output = perf_Bayes_detect(SNR, M, crit)
%% compute Pd and P_FA, PE for Bayes detector
% SNR(SNR<0)=0; 
% crit = th_Bayes;

PD_b=[]; P_FA=[];
for i = 1:length(SNR)
    % right-tail probability=(1-CDF)
    PD_b(i)=1-ncfcdf(crit(i),2,2*M,2*SNR(i));%ncfcdf: non-central CDF of F distr.
    P_FA(i)=1-fcdf(crit(i),2,2*M);%ncfcdf: non-central CDF of F distr.
end

beta = 1-PD_b;
PE = 0.5*(beta + P_FA); % 0.5 is the prior

output = {};
output.Pd= PD_b;
output.PE= PE;
output.PFA= P_FA;

end


function confid = getConfid(SNR, M)
% SNR = signal power/noise power
% SNR --> confid: compute LR (likelihood ratio)
SNR(SNR<0)=0;

L0=[];L1=[];
F_theory = SNR+1; % expected F score
for i = 1:length(SNR)
    L1(i) = ncfpdf(F_theory(i), 2, 2*M, 2*SNR(i));
    L0(i) = fpdf(F_theory(i), 2, 2*M);
end

LR = L1./L0;
confid = LR./(1+LR);
end


function [inter_x, inter_y] = getBayes_threshold(SNR, M)
% SNR = signal power/noise power
% Compute (theory) bayes threshold from PDFs of F Distributions  
SNRdBmax = 20; % 20 dB % maximum value of bayes threshold (intersetion)
SNRmax = 10^(SNRdBmax/10);
x = (0:0.01:SNRmax); % resolution of intersetion

SNR(SNR<0)=0; % ncfpdf accept only SNR (>=0)

N = length(SNR);
inter_x = nan(1,N); inter_y = nan(1,N);
% inter_x = ones(1,N); inter_y = zeros(1,N);
for i = 1:length(SNR)  
    SNR1 = SNR(i);
    if (SNR1~=0) % nonzero
        p1 = ncfpdf(x,2,2*M,2*SNR1);
        p0 = fpdf(x,2,2*M);
        [inter_x(i), inter_y(i)] = find_intersect(p0,p1,x);
    end
end
end


function [inter_x, inter_y] = getBayes_threshold_bf(SNR, M, bf)
% SNR = signal power/noise power
% Compute (theory) bayes threshold from PDFs of F Distributions 
% bf: bayes factor, e.g., bf = 3.
% inter_x (>1): th_BF > 1 (F=SNR+1). any threshold (of F) should > 1.

SNRdBmax = 30; % 20 dB % maximum value of bayes threshold (intersetion)
SNRmax = 10^(SNRdBmax/10);
x = (0:0.01:SNRmax); % resolution of intersetion

SNR(SNR<0)=0; % ncfpdf accept only SNR >=0

N = length(SNR);
inter_x = nan(1,N); inter_y = nan(1,N); % without boundary
% inter_x = ones(1,N); inter_y = zeros(1,N); % set low boundary
for i = 1:length(SNR)  
    SNR1 = SNR(i);
    if (SNR1~=0) % nonzero
        p1 = ncfpdf(x,2,2*M,2*SNR1);
        p0 = fpdf(x,2,2*M);
        [inter_x(i), inter_y(i)] = find_intersect(bf*p0,p1,x);
    end
end
end


function [x, y] = find_intersect(p0,p1,x_axis)
% find intersect of two PDFs
err = abs(p0-p1);
% err2 = flip(err);

dif = diff(err);
% dif2 = diff(err2);

% find where it change from negtive to positive (local minimum)
loc = [];
for i=2:length(dif)
    if dif(i-1)<0 && dif(i)>0
        loc=i;
        break;
    end
end

x = x_axis(loc);
y = 0.5*(p0(loc)+ p1(loc));

if isempty(x) % find no points within searching range
    x = nan;
    y = nan;
end

end


function SNR = SNRdB_Nt2SNR_1Tr(F_dB, Nt)
% from F (dB) of Nt trials to single-trial SNR (signal power/noise power)
SNR = (10.^(F_dB./10)-1)./Nt;

end


function output = get_Pd_SPRT(SNRdB_true, Ns_max, N_MC, M, plt)
% compute Pd of SPRT on each no. of trials

%% (1) generate data
% max sample size of data we observe
% Ns_max = 1000;
% N_MC = 1e3; % 
% % the true parameters, but of course we do not see these values...
% SNRdB_true = -10; % dB

% we see the data generated, dependent on the above values.
SNR_true = 10.^(SNRdB_true./10);
% M = 12; % n of neighbering fins
data = ncfrnd(2, 2*M, 2*SNR_true, N_MC, Ns_max); % H1

% check the Distribution of data
% nbins = 50;
% figure
%     histogram(data, nbins);
%     legend('H_1: F')
%     xlabel('F')
%     ylabel('Count')
%     title('Distribution of the F')

%% (2) compute LLR
step_sig = 10;
Nt = linspace(0,Ns_max,(Ns_max/step_sig)+1);
Nt(1) = []; %starting from 10
    
% compute the max(likelihood) for each round of MC
LLR = [];
parfor i = 1:length(Nt)
% for i = 1:length(Nt)
    data_tmp = data(:, 1:Nt(i));
    p1 = ncfpdf(data_tmp, 2, 2*M, 2*SNR_true);
    p0 = fpdf(data_tmp, 2, 2*M);
    LLR(:,i)= sum(log(p1),2) - sum(log(p0),2);
end

LLRm = mean(LLR,1);
LLRsd = std(LLR,1,1);

%% (3) LLR > Tr (critical value)
lower_lim = -1.56;

upper_lim = log(0.8/0.05); % <<<<<<<<<<<<<<<<<<<<<<<
% upper_lim = log(0.5/0.05);

LLR_tr = (LLR > upper_lim);
LLR_sum = sum(LLR_tr,1); 
Pd = LLR_sum/N_MC;

output = {};
output.Pd = Pd;
output.LLRm = LLRm;
output.LLRsd=LLRsd;
output.Nt = Nt;

%%
if plt==1
    figure
        subplot(211)
        plot(Nt, Pd)
        grid on
        xlabel('# trials')
        ylabel('P_D')
        subplot(212)
        errorbar(Nt, LLRm, LLRsd); hold on
        plot(Nt, upper_lim*ones(size(Nt)), '--k');
        grid on
        xlabel('# trials')
        ylabel('LLR')
end
end



function [LR, Nt, SNR_ML_step] = SPRT_on_F_ML_SNR(data, n, M, SNRdB_lower)
% SNR_true is unknown
% lower bound of single-trial SNR, e.g., SNR_lower = -23 dB
% n=10;% compute every 10 points 
N = length(data);
Nt = linspace(0,N,(N/n)+1);
Nt(1) = [];

LR = []; % log(e) likelihood ratio
SNR_ML_step = [];
for i = 1:length(Nt)
    data_tmp = data(1:Nt(i));
    SNR_ML = ML_snr(data_tmp, SNRdB_lower, M); % power ratio
    p1 = ncfpdf(data_tmp, 2, 2*M, 2*SNR_ML);
    p0 = fpdf(data_tmp, 2, 2*M); 
    LR(i) = sum(log(p1)) - sum(log(p0));
    SNR_ML_step(i)=10*log10(SNR_ML);
end
% LR(1)=nan;
% SNR_ML_step(1)=nan;
end


function SNR_ML = ML_snr(data, SNR_lower, M)
% estimate single-trial SNR(power ratio) by Maximum Likelihood

% (1) define range of SNR of ML
SNR_upper = 20; % dB
%every 1 dB:
SNRdB = linspace(SNR_lower, SNR_upper, SNR_upper-SNR_lower+1); 
% from dB to power ratio:
SNR_x = 10.^(SNRdB./10); % SNR = signal power/noise power

% (2) compute the max(likelihood) for each round of MC
logL = [];
for i = 1:length(SNR_x)
    p1 = ncfpdf(data, 2, 2*M, 2*SNR_x(i));
    logL(i)= sum(log(p1));
end

[~, I] = max(logL);
SNR_ML = SNR_x(I); % power ratio
end





function [LR, Nt] = SPRT_on_Fscore(data, n, M, SNRdB_true)
% SNR_true is known
% n=10;% compute every 10 points 
SNR_true = 10.^(SNRdB_true./10);
N = length(data);
Nt = linspace(0,N,(N/n)+1);
LR = [];
for i = 2:length(Nt)
    data_tmp = data(1:Nt(i));
    p1 = ncfpdf(data_tmp, 2, 2*M, 2*SNR_true);
    p0 = fpdf(data_tmp, 2, 2*M); 
    LR(i) = sum(log(p1)) - sum(log(p0));
end

% compute 1st LR
p1 = ncfpdf(data(1), 2, 2*M, 2*SNR_true);
p0 = fpdf(data(1), 2, 2*M);
LR(1) = sum(log(p1)) - sum(log(p0));

end


function [muHat, sigmaHat, PD, crit_NP] = NP_detector_SNR(SNR, M, logF)
% M = 16; % number of neighboring f bins (noise)
% SNR = 10; % here, SNR is raio (i.e., SNR >0), not dB

%% generate signals
Ntrial = 1e5;
% fix the random number generator
% rng(0,'twister');

mu = 0;
sigma = 1; % standard deviation

x0_k = normrnd(mu,sigma,1,Ntrial) + 1i*normrnd(mu,sigma,1,Ntrial);  % noise
x1_k = sqrt(2*(sigma^2)*SNR) + normrnd(mu,sigma,1,Ntrial) + ...
        1i*normrnd(mu,sigma,1,Ntrial); % signal
    
% complex signals --> complex magnitude (power = magnitude^2) 
x0_k_power = abs(x0_k).^2; 
x1_k_power = abs(x1_k).^2; 

% generate M*Ntrial neighboring bins
x0_k_m = normrnd(mu,sigma,M,Ntrial) + 1i*normrnd(mu,sigma,M,Ntrial); 
x0_k_m_power = abs(x0_k_m).^2; 
Ftest_denominator = mean(x0_k_m_power, 1);

% spectral F-test:
Ftest_score_1 = x1_k_power./Ftest_denominator;% numerator is signal
Ftest_score_0 = x0_k_power./Ftest_denominator;% numerator is noise
   
if logF
    % mean of Ftest_score_1 and 95% CI
    tmp = 10*log10(Ftest_score_1);
    [muHat,sigmaHat,muCI,sigmaCI] = normfit(tmp); 
%     muHat = quantile(tmp,[0.50]);
%     muCI95 = quantile(tmp,[0.05 0.95]);
%     sigmaHat = (muCI95(2)-muCI95(1))/2;
else
    [muHat,sigmaHat,muCI,sigmaCI] = normfit(Ftest_score_1); 
end

% observe values greater than crit_NP only 5% of the time by chance
crit_NP = finv(0.95,2,2*M); %F inverse cumulative distribution function
% crit_NP = finv(0.99,2,2*M); %F inverse cumulative distribution function

PD = sum(Ftest_score_1>crit_NP)/length(Ftest_score_1);
end



function P_D = Pd_F_distribution(P_FA, M, SNR)
% compute theory PD for an F distribution.
% P_FA = 0.05;
% M = 12;
% SNR = signal power/noise power
% from dB to power ratio: SNR = 10.^(SNRdB./10); 

% (example) compute theory PD for Guassian distribution:
% P_FA = 0.05;
% P_D = 0.5*erfc(erfcinv(2*P_FA) - (1/sqrt(2))*sqrt(SNR));

SNR(SNR<0)=0; 
crit_FA = finv(1-P_FA,2,2*M); %F inverse cumulative distribution function
P_D=[];
for i = 1:length(SNR)
    % right-tail probability=(1-CDF)
    P_D(i)=1-ncfcdf(crit_FA,2,2*M,2*SNR(i));%ncfcdf: non-central CDF of F distr.
end
end


function [SNRdB, posteriror_log] = MAP_estimate_SNRofonetrial(N, M, SNRdB_true)
% Maximum A-posterior (MAP) estimation of the SNR of one sigal trial

%% generate data
% sample size of data we observe, trying varying this
% N = 10

% the true parameters, but of course we do not see these values...
% SNRdB_true = 0; % dB
SNR_true = 10.^(SNRdB_true./10);

% we see the data generated, dependent on the above values.
% M = 12; % n of neighbering fins
data = ncfrnd(2, 2*12, 2*SNR_true, 1,N);

% check the Distribution of data
% nbins = 10;
% figure
%     histogram(data, nbins);
%     legend('H_1: F')
%     xlabel('F')
%     ylabel('Count')
%     title('Distribution of the F')

%% Maximum A-posterior (MAP)
% (1) define the prior: uniform in [-30,20]
a = -30; b = 20;
SNRdB = linspace(a,b,501); 
% from dB to power ratio:
SNR_x = 10.^(SNRdB./10); % SNR = signal power/noise power

%Continuous uniform probability density function
p_SNRdB = unifpdf(SNRdB,a,b);  

% (2) compute the likelihood
likelihood_log = [];
for i = 1:length(SNR_x)
    p1 = ncfpdf(data, 2, 2*M, 2*SNR_x(i));
    likelihood_log(i)= sum(log(p1));
end

% (3) compute the posteriror
posteriror_log = likelihood_log + log(p_SNRdB);

end
