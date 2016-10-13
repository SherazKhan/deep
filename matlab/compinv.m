SUBJECTS_DIR=' /cluster/transcend/MRI/WMA/recons';
subj ='hari';
megfile=' /autofs/cluster/transcend/sheraz/data/MEG_EEG/taskforce_1_rest.fif';


Command_fwd=['mne_do_forward_solution --meas /autofs/cluster/transcend/MEG/fix/' subj '/X/' subj '_fix_1_raw.fif --megonly --overwrite ' ...
    ' --spacing  oct-5  --mri ' '/autofs/cluster/transcend/MEG/erm/' subj '/X/' subj '_1-trans.fif    --subject  ' subj  ...
    ' --fwd  /autofs/cluster/transcend/MEG/fix/' subj '/X/' subj '-oct_fix_1-fwd.fif '];

 command = ['mne_do_forward_solution --subject hari  --mindist 5  --src /cluster/transcend/MRI/WMA/recons/hari/bem/hari-volume-7mm-src.fif -- trans /autofs/cluster/transcend/sheraz/data/MEG_EEG/hari-trans.fif  --meas /autofs/cluster/transcend/sheraz/data/MEG_EEG/taskforce_1_rest.fif --bem hari-5120 --megonly --overwrite --fwd meg-vol-7-fwd.fif']
%%
 % 
% 
Command_inv =[' mne_do_inverse_operator --meg --fixed   --noiserank 61 ' ...
    ' --fwd /autofs/eris/p41p3/john/data/MEG_EEG/forward_solution.fif --proj /autofs/cluster/transcend/MEG/fix/' subj '/X/' subj '_fix_1_ecg_proj.fif '...
    '  --proj /autofs/cluster/transcend/MEG/fix/' subj '/X/' subj '_fix_1_eog_proj.fif  '....
    ' --senscov /autofs/cluster/transcend/MEG/fix/' subj '/X/' subj '_erm_1_0.1-144fil-cov.fif   ' ...
   '  '...
   ' --inv /autofs/cluster/transcend/MEG/fix/' subj '/X/' subj '_fix_0.1_144_fil_fixed_new_erm_megreg_0_new_MNE_proj-inv.fif ' ...
  '  '];


%[s,r]=unix(Command_fwd)

%[s,r]=unix(Command_inv)

%a
% 
% command=strcat(' mne_do_forward_solution --subject hari ',['  --mindist 5  --src  '], SUBJECTS_DIR,...
%     '/', SUBJECT, '/bem/volume-7mm-src.fif  --meas  ',megfile, ' --bem sample-5120 --megonly --overwrite --fwd meg-vol-7-fwd.fif');
% 
% 
%  
% 
   [status,result] = unix(Command_inv)
%   
%   fwd=mne_read_forward_solution('meg-vol-7-fwd.fif',1);