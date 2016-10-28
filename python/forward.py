import mne
conductivity = (0.3,)
subject = 'hari'
subjects_dir = '/eris/p41p3/john/data/MEG_EEG'
fname_raw = '/eris/p41p3/john/data/MEG_EEG/meg_only_assr_raw.fif'
trans = '/eris/p41p3/john/data/MEG_EEG/hari-trans.fif'

model = mne.make_bem_model(subject=subject, ico=4,
                           conductivity=conductivity,
                           subjects_dir=subjects_dir)
bem = mne.make_bem_solution(model)
src = mne.setup_source_space(subject, spacing='ico5',
                             subjects_dir=subjects_dir,
                             add_dist=False, overwrite=True)
fwd = mne.make_forward_solution(fname_raw, trans=trans, src=src, bem=bem,fname=None, meg=True, eeg=False, n_jobs=2)

grad_map = mne.sensitivity_map(fwd, ch_type='grad', mode='fixed')
mag_map = mne.sensitivity_map(fwd, ch_type='mag', mode='fixed')

mag_map.plot(hemi='split',subjects_dir=subjects_dir,subject=subject,
                      views=['lateral','medial'], surface='inflated')