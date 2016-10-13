clear all;clc;close all;

%Set default text interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%Add path to the MNE matlab toolbox
addpath('/usr/pubsw/packages/mne/stable/share/matlab/');
addpath('/autofs/eris/p41p3/john/scripts/mne-matlab/matlab');

%Read brainstem response data
assr = fiff_setup_read_raw('/autofs/eris/p41p3/john/data/MEG_EEG/assr_data_john_raw.fif');
[alldata,t] = fiff_read_raw_segment(assr,assr.first_samp,assr.last_samp,1:306);
sensors = alldata(1:306,:);
Fs = assr.info.sfreq;

%Change Nyquist frequency to 300 to make it a multiple of 60Hz so that we
%can use Notch filtering.
newTime = (t(1):1/600:t(end));
data = interp1(t,sensors',newTime,'pchip');

%Use Notch filter to get rid of AC harmonics (60,120,180,...)
d = fdesign.comb('notch','N,BW',10,0.5,600);
Hd = design(d);
%fvtool(Hd)
%Adjust the filter as to make it to zero-phase distortion
y = filtfilt(Hd.numerator,Hd.denominator,data)';
Fs = 600;


%%

ChangeSNR = zeros(102,1);
magindex = (3:3:306);
gradindex = (1:306);
gradindex(magindex) = [];

sensorsProjected=tempProjection_matti(y,[],50,magindex,gradindex);%65 works great for assr


mag_avg_orig = mean(y(magindex,:));
mag_avg_proj = mean(sensorsProjected(magindex,:));

[Y_orig_avg,f]=getFFT(mag_avg_orig,Fs);
Y_proj_avg=getFFT(mag_avg_proj,Fs);
figure;
h(1) = subplot(1,2,1);
plot(f,Y_orig_avg);
ylabel('fT/$\sqrt{Hz}$','rot',0)
xlabel('Unfiltered data')
h(2) = subplot(1,2,2);
plot(f,Y_proj_avg);
xlabel('Filtered data')
%linkaxes(h,'y')
%ylim(h(1),get(h(1),'Ylim'))

%%
% 
% 
% 
% ch = 90;
% [y,f] = getFFT(sensors(ch,:),Fs);
% yp = getFFT(sensorsProjected(ch,:),Fs);
% mydiff = y - yp;
% max(mydiff)
% figure;
% h(1) = subplot(1,2,1); plot(f,y);
% h(2) = subplot(1,2,2); plot(f,yp);
% linkaxes(h,'y')
% ylim(h(1),get(h(1),'Ylim'))
% 
% 
% 
% [Y_orig_avg,f]=getFFT(mag_avg_orig,Fs);
% Y_proj_avg=getFFT(mag_avg_proj,Fs);
% figure; subplot(1,2,1); plot(f,Y_orig_avg); subplot(1,2,2); plot(f,Y_proj_avg)



%Plot original time signal, then plot projected signalfsafs
%figure;plotyy(t,sensors(60,:),t,sensorsProjected(60,:));


k = 100;
h = figure;
subplot(3,6,1)
plot_h2(y(magindex,k));
freezeColors(h)
subplot(3,6,2)
plot_h2(y(magindex,k+1))
subplot(3,6,3)
plot_h2(y(magindex,k+2))
subplot(3,6,4)
plot_h2(y(magindex,k+3))
subplot(3,6,5)
plot_h2(y(magindex,k+4))
subplot(3,6,6)
plot_h2(y(magindex,k+5))
colorbar

subplot(3,6,7)
plot_h2(sensorsProjected(magindex,k))
subplot(3,6,8)
plot_h2(sensorsProjected(magindex,k+1))
subplot(3,6,9)
plot_h2(sensorsProjected(magindex,k+2))
subplot(3,6,10)
plot_h2(sensorsProjected(magindex,k+3))
subplot(3,6,11)
plot_h2(sensorsProjected(magindex,k+4))
subplot(3,6,12)
plot_h2(sensorsProjected(magindex,k+5))

colorbar







