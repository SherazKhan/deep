function [Y,f]=getFFT(y,Fs)


% T = 1/Fs;
% yf = fft(y);
% L = length(y);
% Y = sqrt(1/(L*Fs)*abs(yf(1:ceil(L/2))+1).^2);
% f = 0:Fs/L:Fs/2;


L = length(y);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Ft = fft(y,NFFT);

%Ft = Ft(1:NFFT/2+1);
%Y = Ft/(NFFT*Fs);

 Ft = abs(Ft(1:NFFT/2+1)).^2;
 Y = sqrt(Ft/(NFFT*Fs));

f = Fs/2*linspace(0,1,NFFT/2+1);

%Y=2*abs(Y(1:NFFT/2+1));
% 
% 
% 
% L = length(y);
% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% Y = fft(y,NFFT)/(L*Fs);
% f = Fs/2*linspace(0,1,NFFT/2+1);
% 
% Y=2*abs(Y(1:NFFT/2+1));
