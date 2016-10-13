clear all;clc;
%Add path to the MNE matlab toolbox
addpath('/usr/pubsw/packages/mne/stable/share/matlab/');
%load and read data
%raw=fiff_setup_read_raw('/autofs/space/amiga_001/users/meg/erm1/AC076/1/AC076_erm_1_sss.fif');
raw=fiff_setup_read_raw('/autofs/eris/p41p3/john/MNE-sample-data/MEG/sample/sample_audvis_raw.fif');
infile='/autofs/eris/p41p3/john/MNE-sample-data/MEG/sample/sample_audvis_raw.fif';
[data,times] = fiff_read_raw_segment(raw,raw.first_samp,raw.last_samp,1:306);

ind_mag=3:3:306;
ind_grad=1:306;
ind_grad(ind_mag)=[];

data_grad=data(ind_grad,:);
data_mag=data(ind_mag,:);

%Truncate data to make it feasable for Matlab to handle, take 3000 first
%time measures
Nr = data_grad(:,1:3000);
Np = data_mag(:,1:3000);
%%%

data_mag_filtered = Np - Np*Nr'*inv(Nr*Nr')*Nr;

temporalProjection('/autofs/space/amiga_001/users/meg/erm1/AC076/1/AC076_erm_1_sss.fif',3);
%%
figure
subplot(1,3,1)
plot(times(1:3000),Np(1,:),'r')
hold on;
plot(times(1:3000),data_mag_filtered(1,:))
subplot(1,3,2)
plot(times(1:3000),Np(2,:),'r')
hold on;
plot(times(1:3000),data_mag_filtered(2,:))
subplot(1,3,3)
plot(times(1:3000),Np(3,:),'r')
hold on;
plot(times(1:3000),data_mag_filtered(3,:))
