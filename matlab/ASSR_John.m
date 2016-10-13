clear all;clc;close all;

%Set default text interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%Add path to the MNE matlab toolbox
addpath('/usr/pubsw/packages/mne/stable/share/matlab/');
addpath('/autofs/eris/p41p3/john/scripts/mne-matlab/matlab');
addpath /autofs/eris/p41p3/john/impress/MEG/Matlab_scripts/core_scripts/Private_epochMEG
addpath /autofs/eris/p41p3/john/impress/MEG/Matlab_scripts/TimeFrequencyKosti
addpath /autofs/eris/p41p3/sheraz/matlab_scripts/mne-matlab/matlab

%Read brainstem response data
% dirpath= '/eris/p41p3/john/data/John_ASSR_43Hz_113Hz/'
% fname=[dirpath '43Hz_113Hz_1000HzCarr_raw.fif']
% eventname = [dirpath '43Hz_113Hz_1000HzCarr_saved-raw-eve.fif']

dirpath= '/autofs/eris/p41p3/data/MEG_EEG70/subj_John_02/161005/';
fname1=[dirpath 'ASSR_43_113_197_271_2KHz_raw.fif'];
eventname1 = [dirpath 'ASSR_43_113_197_271_2KHz_raw-eve.fif'];
fname2=[dirpath 'ASSR_43_113_197_271_2KHz_raw-1.fif'];
eventname2 = [dirpath 'ASSR_43_113_197_271_2KHz_raw-1-eve.fif'];
fname3=[dirpath 'ASSR_43_113_197_271_2KHz_raw-2.fif'];
eventname3 = [dirpath 'ASSR_43_113_197_271_2KHz_raw-2-eve.fif'];

% [data_1,time_1] = mne_read_epochs(fname,1,eventname,-0.5,2);
[data_4_1,time_4_1] = mne_read_epochs(fname1,4,eventname1,-0.5,2);
[data_4_2,time_4_2] = mne_read_epochs(fname2,4,eventname2,-0.5,2);
[data_4_3,time_4_3] = mne_read_epochs(fname3,4,eventname3,-0.5,2);

%Append the data into one 3D matrix
l1 = length(data_4_1(1,1,:));
l2 = length(data_4_2(1,1,:));
l3 = length(data_4_3(1,1,:));

data4 = zeros(length(data_4_1(:,1,1)),length(data_4_1(1,:,1)),length(data_4_1(1,1,:)) + length(data_4_2(1,1,:)) + length(data_4_3(1,1,:)));% + size(data_4_2) + size(data_4_3);
data4(:,:,(1:l1)) = data_4_1;
data4(:,:,(l1+1:l1+l2)) = data_4_2;
data4(:,:,(l1+l2+1:end)) = data_4_3;

assrraw = fiff_setup_read_raw([dirpath 't_raw.fif']);
bads = assrraw.info.bads;
badch = [];
for k = 1:length(bads)
    for i = 1:306%length(assrraw.info.ch_names)
        if strcmp(bads{k},assrraw.info.ch_names{i})
            badch = [badch i];
        end
    end
end

%let's just do MEG for now
data = data4(1:306,:,:);


% %[noise,noisetimes] = fiff_read_raw_segment(assrraw,assrraw.first_samp,assrraw.last_samp,1:306);
% assrraw.info.bads = {'MEG1133(275)' 'MEG21430813(85)'};

%Set grad and mag index, remove bad channels.
gradindex = (1:306);
magindex = (3:3:306);
gradindex(magindex) = [];
badindex = badch;
[~,in] = min(abs(gradindex-badindex(1)));
gradindex(in) = [];
[~,in] = min(abs(gradindex-badindex(2)));
gradindex(in) = [];
allindex = (1:306);
allindex(badindex) = [];

Fs = assrraw.info.sfreq;
%mydat = squeeze(data_2(150,(500:1500),:));
meandata = double(mean(data(:,:,:),3));
% % [A,f] = getFFT(meandata(275,:),Fs);
% % plot(f,A)
%Change Nyquist frequency to 300 to make it a multiple of 60Hz so that we
%can use Notch filtering.
% newTime = (time_2(1):1/600:time_2(end));
% data = interp1(time_2,meandata',newTime,'pchip');
% 
% %Use Notch filter to get rid of AC harmonics (60,120,180,...)
% d = fdesign.comb('notch','N,BW',10,0.5,600);
% Hd = design(d);
% %fvtool(Hd)
% %Adjust the filter as to make it to zero-phase distortion
% y = filtfilt(Hd.numerator,Hd.denominator,data)';
% newFs = 600;

% % figure
% % [An,f] = getFFT(y(153,:),newFs);
% % plot(f,An)


%%
y = meandata;
newFs = Fs;

ChangeSNR = zeros(102,1);

sensorsProjected=tempProjection_matti(y,[],90,magindex,gradindex);%65 works great for assr


sig1 = mean(y(magindex,:));
sig2 = mean(sensorsProjected(magindex,:));
sig3 = mean(y(gradindex,:));

sig1 = y(150,:);
sig2 = sensorsProjected(150,:);

[Y_orig_avg,f]=getFFT(sig1,newFs);
Y_proj_avg=getFFT(sig2,newFs);
figure;
h(1) = subplot(1,2,1);
plot(f,Y_orig_avg);
ylabel('fT/$\sqrt{Hz}$','rot',0)
xlabel('Unfiltered data')
h(2) = subplot(1,2,2);
plot(f,Y_proj_avg);
xlabel('Filtered data')
% linkaxes(h,'y')
% ylim(h(1),get(h(1),'Ylim'))

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







