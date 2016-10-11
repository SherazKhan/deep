import mne
from anlffr import spectral
import matplotlib.pylab as plt
plt.ion()

raw_fname = '/eris/p41p3/data/MEG_EEG70/subj_John_02/161005/ASSR_43_113_197_271_2KHz_raw.fif'
raw_bad = '/eris/p41p3/data/MEG_EEG70/subj_John_02/161005/t.fif'

raw = mne.io.read_raw_fif(raw_fname,preload=True)
raw_bad = mne.io.read_raw_fif('/eris/p41p3/data/MEG_EEG70/subj_John_02/161005/t_raw.fif')

raw.info['bads'] = raw_bad.info['bads']
events = mne.find_events(raw, stim_channel='STI101',
                         mask=15, mask_type='and')

event_id = dict(Freq43Hz=1, Freq113Hz=2, Freq197Hz=3, Freq271Hz=4)

tmin, tmax = -0.05, 2

picks = mne.pick_types(raw.info, meg='mag', eeg=False, eog=False, stim=False,
                       exclude='bads')

reject = dict(mag=4e-12)

epochs_mag = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
                    baseline=(-.05, 0), reject=reject, proj=True)

picks = mne.pick_types(raw.info, meg='grad', eeg=False, eog=False, stim=False,
                       exclude='bads')

reject = dict(grad=4000e-13)

epochs_grad = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
                    baseline=(-.05, 0), reject=reject, proj=True)

x = epochs_grad['Freq43Hz'].get_data()
x = x.transpose((1,0,2))
params = dict(Fs=2000,fpass=[5,600],tapers=[2, 3],itc=1)
(plv,f) = spectral.mtcpca(x,params)

plt.plot(f,plv,linewidth = 2)
plt.xlabel('Frequency (Hz)')
plt.grid(True)
plt.show()

x = epochs_mag['Freq113Hz'].get_data()
x = x.transpose((1,0,2))
params = dict(Fs=2000,fpass=[5,600],tapers=[2, 3],itc=1)
(plv,f) = spectral.mtcpca(x,params)

plt.plot(f,plv,linewidth = 2)
plt.xlabel('Frequency (Hz)')
plt.grid(True)
plt.show()

x = epochs_grad['Freq113Hz'].get_data()
x = x.transpose((1,0,2))
params = dict(Fs=2000,fpass=[5,600],tapers=[2, 3],itc=1)
(plv,f) = spectral.mtcpca(x,params)
plt.figure(2)
plt.plot(f,plv,linewidth = 2)
plt.xlabel('Frequency (Hz)')
plt.grid(True)
plt.show()