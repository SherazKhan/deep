clear all;clc;close all;

%Add path to the MNE matlab toolbox
addpath('/usr/pubsw/packages/mne/stable/share/matlab/');
addpath('/autofs/eris/p41p3/john/scripts/mne-matlab/matlab');

%read forward solution
fwd=mne_read_forward_solution('/autofs/eris/p41p3/john/data/MEG_EEG/meg-vol-7-fwd.fif');

%Read real background (resting state) noise
rest = fiff_setup_read_raw('/autofs/eris/p41p3/john/data/MEG_EEG/taskforce_1_rest_filter_raw.fif');


%Find source- and sensorspace (and plot them)
[sensorspace,sourcespace] = FindSensSourceSpace(fwd);

%Gain matrix G
c = struct2cell(fwd);
solstr = c{5};
solcell = struct2cell(solstr);
G = solcell{5};

%Signal. Find grid point in source space closest to zero (approximately the
%location of brainstem)
Fs=250;
dist = sqrt(sum(sourcespace.^2,2));
[vale,ind] = min(dist);
myx = (60:1/Fs:180);
N=length(myx);
sourceinddeep = zeros(3806*3,N);

%Get noise from the rest state data (60-180)
[noise,noisetimes] = fiff_read_raw_segment(rest,rest.first_samp,rest.last_samp,1:306);
%Remove bad channels
% noise(45,:) = zeros(size(noise(45,:)));
% noise(103,:) = zeros(size(noise(103,:)));
% noise(179,:) = zeros(size(noise(179,:)));
% for i = (1:length(rest.info.ch_names))
%     for k = (1:length(rest.info.bads))
%         if rest.info.ch_names{i} == rest.info.bads{k}
%             
%         end
%     end
%     mag_ind = rest.info.ch_names{1}
% end


[~,indstart] = min(abs(noisetimes-60));
[~,indend] = min(abs(noisetimes-180));
noise = noise(:,indstart:indend);
noisetimes = noisetimes(indstart:indend);

%Forward simulation deep brain
DeepFreq = 50;
sourceinddeep((ind-1)*3+1:ind*3,:) = repmat(10^(-8)*sin(2*pi*DeepFreq*myx),3,1);%10 nAm dipole strength
sensinddeep = G*sourceinddeep;

%calculate SNR
totalpowerdeepbrain = sum(sum(sensinddeep.^2))/length(myx);
noisepower = sum(sum(noise.^2))/length(noise);
SNR = totalpowerdeepbrain/noisepower;
dbrms = sqrt(sum(sum(sensinddeep.^2)));
noiserms = sqrt(sum(sum(noise.^2)));

%Superpose signals
sensors = sensinddeep + noise;

ChangeSNR = zeros(306,1);
magindex = (3:3:306);
gradindex = (1:306);
gradindex(magindex) = [];

er = round(rand*102)*3;
sensorsProjected=tempProjection_matti(sensors,[],90);
%plotFFT(sensorsProjected(60,:),250)
%Break point between 1.4 and 1.5 rad!

for chl = (3:3:306)
    
    sig=sensors(chl,:);
    
    
    [Y_orig,f]=getFFT(sig,Fs);
    
    [~,ind10]=min(abs(f-10));
    [~,ind50]=min(abs(f-DeepFreq));
    
    
    
    sig=sensorsProjected(chl,:);
    Y_proj=getFFT(sig,Fs);
    
    SNR_orig=Y_orig(ind50)/Y_orig(ind10);
    SNR_proj=Y_proj(ind50)/Y_proj(ind10);
    %         n = sum(Y_orig)-(Y_orig(ind10)+Y_orig(ind50));
    
    ChangeSNR(chl) = SNR_proj/SNR_orig;
    if ChangeSNR(chl) < 1
        ChangeSNR(chl) = -1/ChangeSNR(chl);
    end
    
    %
    %           figure;plot(f,Y_orig)
    %           figure;plot(f,Y_proj)
    %
    %
    %
    %     sensorsProjected_alog2=temporalProjection_matrix(sensors,30,Fs);
    %
    %     sig=sensorsProjected_alog2(chl,:);
    %     Y_proj_2=getFFT(sig,Fs);
    %    figure;plot(f,Y_proj_2)
    
end


%%

[Y_orig,f]=getFFT(sensors(180,:),Fs);
Y_proj=getFFT(sensorsProjected(180,:),Fs);
breakind = zeros(8,1);
breakind(1) = 1;
[~,breakind(2)]=min(abs(f-5));
[~,breakind(3)]=min(abs(f-15));
[~,breakind(4)]=min(abs(f-25));
[~,breakind(5)]=min(abs(f-35));
[~,breakind(6)]=min(abs(f-45));
[~,breakind(7)]=min(abs(f-55));
breakind(8)=length(Y_orig);
discamppr = zeros(7,1);
discampor = zeros(7,1);
discfr = 0:10:60;
for k = (1:7)
    discamppr(k) = sum(Y_proj(breakind(k):breakind(k+1)));
    discampor(k) = sum(Y_orig(breakind(k):breakind(k+1)));
end
figure;
plotyy(discfr,discamppr,discfr,discampor)
legend projected original

mag_avg_orig = mean(sensors(magindex,:));
mag_avg_proj = mean(sensorsProjected(magindex,:));

[Y_orig_avg,f]=getFFT(mag_avg_orig,Fs);
Y_proj_avg=getFFT(mag_avg_proj,Fs);
figure; subplot(1,2,1); plot(f,Y_orig_avg); subplot(1,2,2); plot(f,Y_proj_avg)


%%
%Plot original time signal, then plot projected signal
figure;plotyy(myx,sensors(120,:),myx,sensorsProjected(120,:));
figure;plotyy(myx,sensinddeep(120,:),myx,sensorsProjected(120,:));
deepcorrs = [];
for i = (3:3:306)
    [deepcorr,p_corr] = corrcoef([sensinddeep(i,:)' sensorsProjected(i,:)']);
    deepcorrs = [deepcorrs; deepcorr(1,2)];
end
legend('Original deep brain','Filtered')

NegDeepcor = find(floor(deepcorrs));
ind_mag_neg = 3*NegDeepcor;

[sensorspace,sourcespace] = FindSensSourceSpace(fwd);

figure;
scatter3(sensorspace(:,1),sensorspace(:,2),sensorspace(:,3),'bo')
hold on
scatter3(sensorspace(NegDeepcor,1),sensorspace(NegDeepcor,2),sensorspace(NegDeepcor,3),'kx')

% figure;plot(myx,10^(-8)*sin(2*pi*50*myx)+noise)
% figure;plot(myx,)

%%
%
%
% P2 = fft(sig);
% P2 = abs(P2/length(sig));
% P1 = P2(1:end/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = 1000*(0:length(myx)/2)/length(myx);
% figure;plot(f,P1)
%
%
%
% figure
%
% P2 = fft(sig);
% P2 = abs(P2/length(sig));
% P1 = P2(1:end/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = 1000*(0:length(myx)/2)/length(myx);
% hold on;plot(f,P1)




%%








% N = raw.last_samp - raw.first_samp + 1;
% %%
% % %signal. Sawtooth shaped, 1000 periods, 20 bootstrap points in each
% % %period
% %
% % tand = (0:1/20:1)';
% % signal = repmat(tand,floor(N/21)-1,1);
% %
% % signal = [signal; zeros(N - length(signal),1)];
% %%
% %Endow cortical signal with pure sinusoid of frequency 10hz, and
% %subcortical with the same at 50hz
% %Sampling frequency - 1000 Hz
% %Length of signal - 166800/1000 = 166.8
% myx = double((0:1:N-1))'/1000;
% signalsurf = sin(2*pi*10*myx);
% signaldeep = sin(2*pi*50*myx);
%
% %%
% %Make signal
% senstimetracesurf = sensindsurf*signalsurf';
% senstimetracedeep = sensinddeep*signaldeep';
%
% signal = senstimetracesurf + senstimetracedeep;
%signal = signalsurf + signaldeep;
sig = signal(75,:);

figure

P2 = fft(sig);
P2 = abs(P2/length(sig));
P1 = P2(1:end/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = 1000*(0:length(myx)/2)/length(myx);
plot(f,P1)
























%%
%Superpose signal with background noise. DON'T DO IT RIGHT NOW

%Filter signal
%Write signal to fiff file, then run it through temporal filter
outfile = 'forward_solution.fif';
infile='/autofs/eris/p41p3/john/MNE-sample-data/MEG/sample/sample_audvis_raw.fif';
raw = fiff_setup_read_raw(infile);
[datasamp,timessamp] = fiff_read_raw_segment(raw,raw.first_samp,raw.last_samp);
[outfid,cals] = fiff_start_writing_raw(outfile,raw.info);
data = signal;
data = [data; datasamp(307:end,:)];
fiff_write_raw_buffer(outfid,data,cals);
fiff_finish_writing_raw(outfid);

temporalProjection('forward_solution.fif',3);

raw=fiff_setup_read_raw('forward_solution-proj.fif');
[data,times] = fiff_read_raw_segment(raw,raw.first_samp,raw.last_samp,1:306);

ind_mag=3:3:306;
ind_grad=1:306;
ind_grad(ind_mag)=[];

data_grad=data(ind_grad,:);
data_mag=data(ind_mag,:);




%%
sig2 = data(152,:);
figure;

P2 = fft(sig2);
P2 = abs(P2/length(sig2));
P1 = P2(1:end/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = 1000*(0:length(myx)/2)/length(myx);
plot(f,P1)



%Do inverse solution



