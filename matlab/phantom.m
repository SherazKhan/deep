clear all; clc; close all;

addpath('/usr/pubsw/packages/mne/stable/share/matlab/')

indMag = (3:3:306);

% rawmat = [];
% for i = (1:32)
%     filename = ['/autofs/eris/p41p3/john/data/phantom/' num2str(i) '.fif'];
%     rawmat = [rawmat; fiff_setup_read_raw(filename)];
% end
% 
% %Deep brain dipole (db); dipole nr 28, located at (0,31.5,12.7), r=34.0
% %Cortical dipole (cr); dipole nr 29, located at (0,13.9,62.4), r=64.0
% 
% dbraw = rawmat(28);
% crraw = rawmat(29);

%Processed data
dbraw = fiff_setup_read_raw('/autofs/eris/p41p3/john/data/phantom/procdata28_raw.fif');
crraw = fiff_setup_read_raw('/autofs/eris/p41p3/john/data/phantom/procdata29_raw.fif');

[dbdata,dbtimes] = fiff_read_raw_segment(dbraw,dbraw.first_samp,dbraw.last_samp,1:306);
[crdata,crtimes] = fiff_read_raw_segment(crraw,crraw.first_samp,crraw.last_samp,1:306);

% %to separate them and evaluate the efficiency of the algorithm, we can
% %phase shift them by 90 degrees, making them orthogonal in a Hilbert space
% %(and apparently linearly independent), T is approximately 0.05, so 1000
% %iterations should be more than enough. Choose channel 68; MEG0
% scprod = zeros(200,1);
% for k = (1:200)
%     scprod(k) = dbdata(68,:)*crdata(68,:)';
%     dbdata = [dbdata(:,end) dbdata(:,1:end-1)];
%     k
% end


%The signals have different average value, so translate them so that the
%average value is zero for each channel in both db and cr. important also
%for PCA (when channels are variables and time instants seen as samples.)
for i = (1:306)
    dbdata(i,:) = dbdata(i,:) - mean(dbdata(i,:));
    crdata(i,:) = crdata(i,:) - mean(crdata(i,:));
end



%translate db signal by 100 time units
dbdata = [dbdata(:,end-100:end) dbdata(:,1:end-101)];

% figure
% plot(dbtimes,dbdata(68,:));
% figure
% plot(crtimes,crdata(68,:));
% %plot combined signal
combsign = dbdata + crdata;
% figure
% plot(dbtimes,combsign(68,:));


%The dipoles are triggered with two 0.05s periods (so 0.1s of activation)
%and then a pause of 0.25s, so the period taken together is 0.35s. db
%starts at 16.26, so periods of signal are 16.26-16.36; 16.61-16.71.
%Cortical starts at 16.13, so 16.13-16.23; 16.48-16.58

%To evaluate degree of SNR increase, get signal energy for 80 first
%periods
[~,indCortStart] = min(abs(crtimes-16.13));
indCortPeriod = (indCortStart:indCortStart+100);
indCortSign = [indCortPeriod];
for i = 1:80
    indCortSign = [indCortSign indCortPeriod+i*350]; 
end

[~,indDbStart] = min(abs(dbtimes-16.26));
indDbPeriod = (indDbStart:indDbStart+100);
indDbSign = [indDbPeriod];
for i = 1:80
    indDbSign = [indDbSign indDbPeriod+i*350]; 
end

%Energy of all channels
deepEnergy = sum(sum(combsign(indMag,indDbSign).^2));
cortEnergy = sum(sum(combsign(indMag,indCortSign).^2));



%Filter
dataProj=tempProjection_matti_old(combsign,[],0.2);

%Check new signal energy
deepEnergyProj = sum(sum(dataProj(indMag,indDbSign).^2));
cortEnergyProj = sum(sum(dataProj(indMag,indCortSign).^2));

SNR_orig = deepEnergy/cortEnergy
SNR_proj = deepEnergyProj/cortEnergyProj
SNR_inc = SNR_proj/SNR_orig

% %Check that deepbrain signal hasn't been canceled
chl = 90;
figure
plot((1:length(indDbSign)),combsign(chl,indDbSign)); hold on
plot((1:length(indDbSign)),dataProj(chl,indDbSign),'--');
legend orig proj
[deep_corr,p_corr] = corrcoef([combsign(chl,indDbSign)' dataProj(chl,indDbSign)'])
title('Deep')

figure
plot((1:length(indDbSign)),combsign(chl,indCortSign)); hold on
plot((1:length(indDbSign)),dataProj(chl,indCortSign),'--');
legend orig proj
title('Cort')

deepcorrs = [];
for i = (3:3:306)
    [deepcorr,p_corr] = corrcoef([combsign(i,indDbSign)' dataProj(i,indDbSign)']);
    deepcorrs = [deepcorrs; deepcorr(1,2)];
end

figure; hist(deepcorrs);

%[cortical_corr,p_corr] = corrcoef([combsign(chl,indCortSign)' dataProj(chl,indCortSign)'])
%[deep_corr,p_deep] = corrcoef([combsign(chl,indDbSign)' dataProj(chl,indDbSign)'])

% figure
% dbtot = sum(dbdata);
% crtot = sum(crdata);
% plot(dbtimes,dbdata(69,:));
% figure
% plot(crtimes,crdata(69,:),'y');
% 
% figure
% plot(dbtimes,combsign(90,:));
% 

% mixdataFiltered = tempProjection_matti(mixdata);
% dbdataFiltered = tempProjection_matti(dbdata);
% crdataFiltered = tempProjection_matti(crdata);




