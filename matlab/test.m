

clear all
raw=fiff_setup_read_raw('/autofs/eris/p41p3/john/MNE-sample-data/MEG/sample/sample_audvis_raw.fif');
[data,times] = fiff_read_raw_segment(raw,raw.first_samp,raw.last_samp);

[outfid,cals] = fiff_start_writing_raw('test.fif',raw.info);
fiff_write_raw_buffer(outfid,data,cals);
fiff_finish_writing_raw(outfid);

temporalProjection('test.fif',3);

%%

outfile = 'forward_solution.fif';
infile='/autofs/eris/p41p3/john/MNE-sample-data/MEG/sample/sample_audvis_raw.fif';
raw2 = fiff_setup_read_raw(infile);
[datasamp,timessamp] = fiff_read_raw_segment(raw2,raw2.first_samp,raw2.last_samp);
[outfid2,cals2] = fiff_start_writing_raw(outfile,raw2.info);
data2 = zeros(size(data));

fiff_write_raw_buffer(outfid2,data2,cals2);
fiff_finish_writing_raw(outfid);

%%


