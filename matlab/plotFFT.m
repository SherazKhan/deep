function plotFFT(y,Fs)




L = length(y);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Ft = fft(y,NFFT);
Ft = abs(Ft(1:NFFT/2+1)).^2;
Y = sqrt(Ft/(L*Fs));
f = Fs/2*linspace(0,1,NFFT/2+1);


figure;plot(f,Y)
