function SNR = get_spectrum_F(f, Y_f, nbin, r)
% faster method for computing F of indeger frequency [1,2,3 ... N-1].
% each target intiger f is aligned to a f bin; 
% The noise is selected as +/- r Hz around a target intiger f, exclude the
% target f bin. (r = 0.5 Hz).
% nbin = 12; % f resolusion = 1/12 Hz

f_set = 1: (max(f)-1); % the last f cannot be computed without neighers
SNR =[];
for i=1:length(f_set)  
    fn = f_set(i);
    k = nbin*fn +1; 
    % for checking
%     disp (f(k));
    mY = (Y_f(k)); 
    nr = floor(r* nbin); % 0.5*12 = 6 (f bins)
    range = Y_f([[k-nr : k-1], [k+1: k+nr]]);
    mY2 = mean(range);
    SNR(i) = (mY/mY2)^2; % uV^2
    % 10*log(mY/mY2) is not clear to show peaks 
%     SNR(i) = 10*log(mY/mY2); % dB. Note: log(1) = 0. 
end



