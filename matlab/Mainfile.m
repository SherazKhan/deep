clear all
%Add path to the MNE matlab toolbox
addpath('/usr/pubsw/packages/mne/stable/share/matlab/');
%load and read data
raw=fiff_setup_read_raw('/autofs/space/amiga_001/users/meg/erm1/AC076/1/AC076_erm_1_sss.fif');
[data2,times] = fiff_read_raw_segment(raw,raw.first_samp,raw.last_samp,1:306);

temporalProjection('/autofs/space/amiga_001/users/meg/erm1/AC076/1/AC076_erm_1_sss.fif');
load('tempprojdata');
%%
ind_mag=3:3:306;
ind_grad=1:306;
ind_grad(ind_mag)=[];

data_grad=data2(ind_grad,:);
data_mag=data2(ind_mag,:);

%Truncate data to make it feasable for Matlab to handle, take 3000 first
%time measures
Nr = data_grad(:,1:3000);
Np = data_mag(:,1:3000);
%%%

data_mag_filtered = Np - Np*Nr'*inv(Nr*Nr')*Nr;

%temporalProjection('/autofs/space/amiga_001/users/meg/erm1/AC076/1/AC076_erm_1_sss.fif',3);

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
subplot(1,3,1)
plot(times(1:3000),data(3,:),'k:')
subplot(1,3,2)
plot(times(1:3000),data(6,:),'k:')
subplot(1,3,3)
plot(times(1:3000),data(9,:),'k:')
%fiff_finish_writing_raw(outfid);

corrcoef(data(3,:)',data_mag_filtered(1,:)')
corrcoef(data(6,:)',data_mag_filtered(2,:)')
corrcoef(data(9,:)',data_mag_filtered(3,:)')

corrcoef(data(3,:)',Np(1,:)')
corrcoef(data(6,:)',Np(2,:)')
corrcoef(data(9,:)',Np(3,:)')
