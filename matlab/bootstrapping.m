clear all;clc;close all;

%Add path to the MNE matlab toolbox
addpath('/usr/pubsw/packages/mne/stable/share/matlab/');
addpath('/autofs/eris/p41p3/john/scripts/mne-matlab/matlab');

%read forward solution
fwd=mne_read_forward_solution('/autofs/eris/p41p3/john/data/MEG_EEG/meg-vol-7-fwd.fif');

%Read real background (resting state) noise
rest = fiff_setup_read_raw('/autofs/eris/p41p3/john/data/MEG_EEG/taskforce_1_rest_filter_raw.fif');

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
myx = (60:1/Fs:300);
N=length(myx);
sourceinddeep = zeros(3806*3,N);
magindex = (3:3:306);
gradindex = (1:306);
gradindex(magindex) = [];

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
[~,indend] = min(abs(noisetimes-300));
noise = noise(:,indstart:indend);
noisetimes = noisetimes(indstart:indend);

%Find deep source locations
point2 = [0 0 0.025];
[vale,ind2] = min(sqrt(sum((repmat(point2,length(sourcespace),1)-sourcespace)'.^2)));
point3 = [0 0 0.05];
[vale,ind3] = min(sqrt(sum((repmat(point3,length(sourcespace),1)-sourcespace)'.^2)));
point4 = [0 0 0.075];
[vale,ind4] = min(sqrt(sum((repmat(point4,length(sourcespace),1)-sourcespace)'.^2)));
point5 = [0 0 0.10];
[vale,ind5] = min(sqrt(sum((repmat(point5,length(sourcespace),1)-sourcespace)'.^2)));
[valesurf,indsurf] = max(sourcespace(:,3));


%Forward simulation deep brain
DeepFreq1 = 60;
sourceinddeep((ind-1)*3+1:ind*3,:) = repmat(10^(-8)*sin(2*pi*DeepFreq1*myx),3,1);%10 nAm dipole strength
% DeepFreq2 = 70;
% sourceinddeep((ind2-1)*3+1:ind2*3,:) = repmat(10^(-8)*sin(2*pi*DeepFreq2*myx),3,1);%10 nAm dipole strength
% DeepFreq3 = 80;
% sourceinddeep((ind3-1)*3+1:ind3*3,:) = repmat(10^(-8)*sin(2*pi*DeepFreq3*myx),3,1);%10 nAm dipole strength
% DeepFreq4 = 90;
% sourceinddeep((ind4-1)*3+1:ind4*3,:) = repmat(10^(-8)*sin(2*pi*DeepFreq4*myx),3,1);%10 nAm dipole strength
% DeepFreq5 = 100;
% sourceinddeep((ind5-1)*3+1:ind5*3,:) = repmat(10^(-8)*sin(2*pi*DeepFreq5*myx),3,1);%10 nAm dipole strength
% DeepFreq6 = 110;
% sourceinddeep((indsurf-1)*3+1:indsurf*3,:) = repmat(10^(-8)*sin(2*pi*DeepFreq6*myx),3,1);%10 nAm dipole strength

sensinddeep = G*sourceinddeep;
N = 100;%Number of bootstrapping samples
T = floor(length(sourceinddeep)/2);
X = floor(T*rand(N,1));
SNR_orig_mat = [];
SNR_proj_mat = [];
SNRchange = [];
n = 1;

for k = (1:N)
    Tstart = floor(180*rand)+60;
    Tend = Tstart + 60;
    [~,indstart] = min(abs(myx-Tstart));
    [~,indend] = min(abs(myx-Tend));
    
    %Superpose signals
    noisetmp = noise(:,indstart:indend);
    sensinddeeptmp = sensinddeep(:,indstart:indend);
    sensors = noisetmp + sensinddeeptmp; 
    
    %calculate original SNR
    SNR_orig = [];
    for i = (1:length(magindex))
        [origFFT,f] = getFFT(sensors(magindex(i),:),Fs);
        [~,ind1] = min(abs(50-f));
        [~,ind55] = min(abs(55-f));
        source1 = origFFT(ind1);
        noiseMedian = median(origFFT(1:ind55));
        SNR_orig = [SNR_orig source1/noiseMedian];
    end
    SNR_orig_mat = [SNR_orig_mat median(SNR_orig)];
    
    data_filtered=tempProjection_matti(sensors,[],90,magindex,gradindex);%Angle;1.45
    
    %calculate filtered SNR
    SNR_proj = [];
    for i = (1:length(magindex))
        [projFFT,f] = getFFT(data_filtered(magindex(i),:),Fs);
        [~,ind1] = min(abs(50-f));
        [~,ind55] = min(abs(55-f));
        source1 = projFFT(ind1);
        noiseMedian = median(projFFT(1:ind55));
        SNR_proj = [SNR_proj source1/noiseMedian];
    end
    SNR_proj_mat = [SNR_proj_mat median(SNR_proj)];
% 
% mag_avg_orig = mean(sensors(magindex,:));
% mag_avg_proj = mean(data_filtered(magindex,:));
% 
%     [Y_orig_avg,f]=getFFT(mag_avg_orig,Fs);
%     [~,ind50] = min(abs(f-50));
%     SNRor = Y_orig_avg(ind50)/mean(Y_orig_avg);
%     Y_proj_avg=getFFT(mag_avg_proj,Fs);
%     [~,ind50] = min(abs(f-50));
%     SNRpr = Y_proj_avg(ind50)/mean(Y_proj_avg);
%    SNRchange = [SNRchange SNRpr/SNRor];


    disp([num2str(n/N*100) '%'])
    n = n+1;
end

SNRchange = SNR_proj_mat./SNR_orig_mat
%Get error estimate for expected value of SNR increase and SNR distribution
histogram(SNRchange,20)
mu = mean(SNRchange)
var = mean(SNRchange.^2)-mean(SNRchange)^2
mu_errorlimit = sqrt(var/length(SNRchange))

%Plot helmet (snapshots of spatial distribution of source space), one
%deepbrain cycle
k = 100;
h = figure;
subplot(3,6,1)
plot_h2(sensors(magindex,k));
freezeColors(h)
subplot(3,6,2)
plot_h2(sensors(magindex,k+1))
subplot(3,6,3)
plot_h2(sensors(magindex,k+2))
subplot(3,6,4)
plot_h2(sensors(magindex,k+3))
subplot(3,6,5)
plot_h2(sensors(magindex,k+4))
subplot(3,6,6)
plot_h2(sensors(magindex,k+5))
colorbar

subplot(3,6,7)
plot_h2(data_filtered(magindex,k))
subplot(3,6,8)
plot_h2(data_filtered(magindex,k+1))
subplot(3,6,9)
plot_h2(data_filtered(magindex,k+2))
subplot(3,6,10)
plot_h2(data_filtered(magindex,k+3))
subplot(3,6,11)
plot_h2(data_filtered(magindex,k+4))
subplot(3,6,12)
plot_h2(data_filtered(magindex,k+5))

colorbar

subplot(3,6,13)
plot_h2(sensinddeep(magindex,k))
subplot(3,6,14)
plot_h2(sensinddeep(magindex,k+1))
subplot(3,6,15)
plot_h2(sensinddeep(magindex,k+2))
subplot(3,6,16)
plot_h2(sensinddeep(magindex,k+3))
subplot(3,6,17)
plot_h2(sensinddeep(magindex,k+4))
subplot(3,6,18)
plot_h2(sensinddeep(magindex,k+5))

colorbar


%%
sensorsProjected = data_filtered;


% [Y_orig,f]=getFFT(sensors(180,:),Fs);
% Y_proj=getFFT(sensorsProjected(180,:),Fs);
% breakind = zeros(8,1);
% breakind(1) = 1;
% [~,breakind(2)]=min(abs(f-5));
% [~,breakind(3)]=min(abs(f-15));
% [~,breakind(4)]=min(abs(f-25));
% [~,breakind(5)]=min(abs(f-35));
% [~,breakind(6)]=min(abs(f-45));
% [~,breakind(7)]=min(abs(f-55));
% breakind(8)=length(Y_orig);
% discamppr = zeros(7,1);
% discampor = zeros(7,1);
% discfr = 0:10:60;
% for k = (1:7)
%     discamppr(k) = sum(Y_proj(breakind(k):breakind(k+1)));
%     discampor(k) = sum(Y_orig(breakind(k):breakind(k+1)));
% end
% figure;
% plotyy(discfr,discamppr,discfr,discampor)
% legend projected original

mag_avg_orig = mean(sensors(magindex,:));
mag_avg_proj = mean(sensorsProjected(magindex,:));

[Y_orig_avg,f]=getFFT(mag_avg_orig,Fs);
[~,ind50] = min(abs(f-50));
SNRor = Y_orig_avg(ind50)/mean(Y_orig_avg);
Y_proj_avg=getFFT(mag_avg_proj,Fs);
[~,ind50] = min(abs(f-50));
SNRpr = Y_proj_avg(ind50)/mean(Y_proj_avg);
figure; subplot(1,2,1); plot(f,Y_orig_avg); subplot(1,2,2); plot(f,Y_proj_avg)
figure;
h(1) = subplot(1,2,1); plot(f,Y_orig_avg);
h(2) = subplot(1,2,2); plot(f,Y_proj_avg);
SNRrate = SNRpr/SNRor

linkaxes(h,'y')
ylim(h(1),get(h(1),'Ylim'))


% ch = 75;
% [y,f] = getFFT(sensors(ch,:),Fs);
% yp = getFFT(sensorsProjected(ch,:),Fs);
% figure;
% h(1) = subplot(1,2,1); plot(f,y);
% h(2) = subplot(1,2,2); plot(f,yp);
% linkaxes(h,'y')
% ylim(h(1),get(h(1),'Ylim'))


%%
%Plot original time signal, then plot projected signal
myx = myx(x:x+T);
figure;plotyy(myx,sensors(120,:),myx,sensorsProjected(120,:));
figure;plotyy(myx,sensinddeep(75,:),myx,sensorsProjected(75,:));
deepcorrs = [];
for i = (3:3:306)
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



