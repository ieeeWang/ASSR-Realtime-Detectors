function [P1, f] = get_fftAmpSpec(X, fs) 
% get the Amplitude Spectrum of X
% f resolusion =  1/t, t is duration (sec) of X
% OUTPUT
% P1 - |Y(f)|
% P2 - |Y(f)|^2


Y = fft(X); % same with Y = fft(X, lenth(X)); 
L = length(X);
P_abs = abs(Y/L);
P1 = P_abs(1 : floor(L/2+1)); % first half including the median point
P1(2:end-1) = 2*P1(2:end-1);
f = fs/2*(0:(2/L):1); % same with f = Fs*(0:(L/2))/L;
% f resolusion = fs/L = 1/t
P2 = P1.^2; 

% figure
%     plot(f, P1) 
%     title('Single-Sided Amplitude Spectrum of X(t)')
%     xlabel('f (Hz)')
%     ylabel('|P1(f)|')