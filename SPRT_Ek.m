clc; clear; close all;
% compute expected stop time of SPRT
%%
SNRdB = linspace(-30,20,51); 
% from dB to power ratio:
SNR = 10.^(SNRdB./10); % SNR = signal power/noise power

%%  
SNRdBmax = 40; % dB % maximum value of bayes threshold (intersetion)
SNRmax = 10^(SNRdBmax/10);
x = (0:0.01:SNRmax); % resolution of intersetion

N = length(SNR);
M = 12;
tic % Elapsed time is 119 s
KL0 = []; KL1 = []; 
for i = 1:length(SNR)  
    p1 = ncfpdf(x,2,2*M,2*SNR(i));
    p0 = fpdf(x,2,2*M);
    KL0(i) = sum(p0 .* (log(p0)-log(p1)));
    KL1(i) = sum(p1 .* (log(p1)-log(p0)));
end
toc

%%
% Function Handle , % see formulas in appendix
E0_up = @(a,b) (1-a)*log((1-a)/b) - a*log((1-b)/a);
E1_up = @(a,b) (1-b)*log((1-b)/a) - b*log((1-a)/b);
 
a = 0.05; b = 0.2;
figure % KL divergence
subplot(211)
    plot(SNRdB, KL0); hold on
    plot(SNRdB, KL1,'--'); 
    legend ({'KL_0','KL_1'},'Location','northeast')
subplot(212)
%     plot(SNRdB, E0_up./KL0); hold on
%     plot(SNRdB, E1_up./KL1,'--');      
    semilogy(SNRdB, E0_up(a,b)./KL0); hold on
    semilogy(SNRdB, E1_up(a,b)./KL1,'--'); 
    ylim([1, 10^3])
    grid on   
    legend ({'E_0','E_1'},'Location','northeast')
    xlabel('SNR (dB)')
    ylabel('Steps')

% plot stop time of SPRT
figure    
    semilogy(SNRdB, E1_up(a,0.1)./KL1,'-*'); hold on
    semilogy(SNRdB, E1_up(a,0.2)./KL1,'-o');
    semilogy(SNRdB, E1_up(a,0.5)./KL1,'-d');
    ylim([1, 10^3])
    grid on   
    legend ({'P_D = 0.9','P_D = 0.8','P_D = 0.5'},'Location','northeast')
    xlabel('SNR (dB)')
    ylabel('E(K)')
    
    
    