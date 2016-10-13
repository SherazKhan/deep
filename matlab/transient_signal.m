clear all;clc;close all;

%Set default text interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

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

%Set grad and mag index, remove bad channels.
gradindex = (1:306);
magindex = (3:3:306);
gradindex(magindex) = [];
badindex = [45 103 179];
magindex(45/3) = [];
[~,in] = min(abs(gradindex-badindex(2)));
gradindex(in) = [];
[~,in] = min(abs(gradindex-badindex(3)));
gradindex(in) = [];
allindex = (1:306);
allindex(badindex) = [];


%Get noise from the rest state data (60-180)
allchannels = fiff_read_raw_segment(rest,rest.first_samp,rest.last_samp);
[noise,noisetimes] = fiff_read_raw_segment(rest,rest.first_samp,rest.last_samp,1:306);


%Signal. Find grid point in source space closest to zero (approximately the
%location of brainstem)
Fs=250;
dist = sqrt(sum(sourcespace.^2,2));
[vale,ind] = min(dist);
myx = noisetimes;%(60:1/Fs:180);
N=length(myx);
sourceinddeep = zeros(3806*3,N);


% [~,indstart] = min(abs(noisetimes-60));
% [~,indend] = min(abs(noisetimes-180));
% noise = noise(:,indstart:indend);
% noisetimes = noisetimes(indstart:indend);

%Find deep source locations
point2 = [0 0 0.025];
[vale,ind2] = min(sqrt(sum((repmat(point2,length(sourcespace),1)-sourcespace)'.^2)));
point3 = [0 0 0.05];
[vale,ind3] = min(sqrt(sum((repmat(point3,length(sourcespace),1)-sourcespace)'.^2)));
point4 = [0 0 0.075];
[vale,ind4] = min(sqrt(sum((repmat(point4,length(sourcespace),1)-sourcespace)'.^2)));
point5 = [0 0 0.10];
[vale,ind5] = min(sqrt(sum((repmat(point5,length(sourcespace),1)-sourcespace)'.^2)));

%Forward simulation deep brain - put it on/off at random times to get
%transient signals.
DeepFreq1 = 10;
T = 2*DeepFreq1^-1;
N = floor((noisetimes(end)-noisetimes(1))/(8*DeepFreq1^-1));
Dt = (0:1/Fs:T*4-1/Fs);
R = [];
for k = (1:N)
    Ts = 3*T*rand;
    Tend = Ts + T;
    sig = heaviside(Dt-Ts) - heaviside(Dt-Tend);
    R = [R sig];
end
R = [R ones(1,length(noisetimes)-length(R))];
%%

sourceinddeep((ind-1)*3+1:ind*3,:) = repmat(10^(-8)*R.*sin(2*pi*DeepFreq1*myx),3,1);%10 nAm dipole strength
DeepFreq2 = 20;
sourceinddeep((ind2-1)*3+1:ind2*3,:) = repmat(10^(-8)*R.*sin(2*pi*DeepFreq2*myx),3,1);%10 nAm dipole strength
DeepFreq3 = 30;
sourceinddeep((ind3-1)*3+1:ind3*3,:) = repmat(10^(-8)*R.*sin(2*pi*DeepFreq3*myx),3,1);%10 nAm dipole strength
DeepFreq4 = 40;
sourceinddeep((ind4-1)*3+1:ind4*3,:) = repmat(10^(-8)*R.*sin(2*pi*DeepFreq4*myx),3,1);%10 nAm dipole strength
DeepFreq5 = 50;
sourceinddeep((ind5-1)*3+1:ind5*3,:) = repmat(10^(-8)*R.*sin(2*pi*DeepFreq5*myx),3,1);%10 nAm dipole strength

sensinddeep = G*sourceinddeep;
%sensinddeep([45 103 147],:) = zeros(size(sensinddeep([45 103 147],:)));

% %calculate SNR
% totalpowerdeepbrain = sum(sum(sensinddeep.^2))/length(myx);
% noisepower = sum(sum(noise.^2))/length(noise);
% SNR = totalpowerdeepbrain/noisepower;
% dbrms = sqrt(sum(sum(sensinddeep.^2)));
% noiserms = sqrt(sum(sum(noise.^2)));

%Superpose signals
sensors = sensinddeep + noise;

% ChangeSNR = zeros(306,1);
% magindex = (3:3:306);
% gradindex = (1:306);
% gradindex(magindex) = [];

%er = round(rand*102)*3;

sensorsProjected=tempProjection_matti(sensors,[],90,magindex,gradindex);%Angle;90

%plotFFT(sensorsProjected(60,:),250)
%Break point between 1.4 and 1.5 rad!

% for chl = (3:3:306)
%     
%     sig=sensors(chl,:);
%     
%     
%     [Y_orig,f]=getFFT(sig,Fs);
%     
%     [~,ind10]=min(abs(f-10));
%     [~,ind50]=min(abs(f-DeepFreq));
%     
%     
%     
%     sig=sensorsProjected(chl,:);
%     Y_proj=getFFT(sig,Fs);
%     
%     SNR_orig=Y_orig(ind50)/Y_orig(ind10);
%     SNR_proj=Y_proj(ind50)/Y_proj(ind10);
%     %         n = sum(Y_orig)-(Y_orig(ind10)+Y_orig(ind50));
%     
%     ChangeSNR(chl) = SNR_proj/SNR_orig;
%     if ChangeSNR(chl) < 1
%         ChangeSNR(chl) = -1/ChangeSNR(chl);
%     end
%     
%     %
%     %           figure;plot(f,Y_orig)
%     %           figure;plot(f,Y_proj)
%     %
%     %
%     %
%     %     sensorsProjected_alog2=temporalProjection_matrix(sensors,30,Fs);
%     %
%     %     sig=sensorsProjected_alog2(chl,:);
%     %     Y_proj_2=getFFT(sig,Fs);
%     %    figure;plot(f,Y_proj_2)
%     
% end


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
figure;
ground_truth = getFFT(mean(sensinddeep(magindex,:)),Fs);
h(1) = subplot(1,3,1); plot(f,10^15*Y_orig_avg);
ylabel('fT/$\sqrt{Hz}$','rot',0)
xlabel('Unfiltered data')
h(2) = subplot(1,3,2); plot(f,10^15*Y_proj_avg);
xlabel('Filtered data')
h(3) = subplot(1,3,3); plot(f,10^15*ground_truth);
xlabel('Pure deep brain signal')
linkaxes(h,'y')
set(h(1),'Ylim',[0 100])
set(h,'Xlim',[0 100])

ch = 90;
[y,f] = getFFT(sensors(ch,:),Fs);
yp = getFFT(sensorsProjected(ch,:),Fs);
mydiff = y - yp;
max(mydiff)
figure;
h(1) = subplot(1,2,1); plot(f,y);
h(2) = subplot(1,2,2); plot(f,yp);
linkaxes(h,'y')
ylim(h(1),get(h(1),'Ylim'))


%%
%Plot original time signal, then plot projected signal
figure;plotyy(myx,sensors(120,:),myx,sensorsProjected(120,:));
figure;plotyy(myx,sensinddeep(75,:),myx,sensorsProjected(75,:));
deepcorrs = [];
for i = magindex
    [deepcorr,p_corr] = corrcoef([sensinddeep(i,:)' sensorsProjected(i,:)']);
    deepcorrs = [deepcorrs; deepcorr(1,2)];
end
legend('Original deep brain','Filtered')

NegDeepcor = find(floor(deepcorrs));
ind_mag_neg = 3*NegDeepcor;
figure; hist(deepcorrs);

[sensorspace,sourcespace] = FindSensSourceSpace(fwd);

figure;
scatter3(sensorspace(:,1),sensorspace(:,2),sensorspace(:,3),'bo')
hold on
scatter3(sensorspace(NegDeepcor,1),sensorspace(NegDeepcor,2),sensorspace(NegDeepcor,3),'kx')


% %%
% %Write to fif file

%Filter signal
%Write signal to fiff file, then run it through temporal filter

%Remove bad channels from the output file
rest.info.nchan = 306;
rest.info.chs(badindex) = [];
rest.info.ch_names(badindex) = [];
rest.cals(badindex) = [];
sensors = sensors(allindex,:);
sensorsProjected = sensorsProjected(allindex,:);

%Print
outfile = 'filtered_data.fif';
[outfid,cals] = fiff_start_writing_raw(outfile,rest.info);
data = [sensorsProjected; allchannels(end-2:end,:)];
fiff_write_raw_buffer(outfid,data,cals);
fiff_finish_writing_raw(outfid);

% 
% 
% 
% %Do inverse solution



