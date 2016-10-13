function data=temporalProjection(infile,window_length,proj_chs,outfile)

% infile = input fif file
% outfile= output fif file
% window_length = window length in time (can be full length of the file)

if ~exist('outfile','var')
    [PATHSTR,NAME,EXT] = fileparts(infile);
    if ~isempty(PATHSTR)
        outfile=[PATHSTR '/' NAME '-proj' EXT];
    else
        outfile=[ NAME '-proj' EXT];
    end
end

if ~exist('proj_chs','var')
  proj_chs = 1:306;
  proj_chs(3:3:306)=[];
end

raw = fiff_setup_read_raw(infile);
nsamp = raw.last_samp - raw.first_samp;
data_length = nsamp/raw.info.sfreq;

if ~exist('window_length','var')
    window_length=nsamp/raw.info.sfreq;
end


nwin = round(data_length/window_length);  % Number of time windows
                                          % for the temporal projection
nsamp_win = nsamp/nwin;  % Number of samples in one window
ch_kinds=zeros(length(raw.info.ch_names),1);
for j = 1:raw.info.nchan
  ch_kinds(j) = raw.info.chs(j).kind;
end
meg_chs = find(ch_kinds==1);

nProj = length(proj_chs);

[outfid,cals] = fiff_start_writing_raw(outfile,raw.info);
for j = 1:nwin
  from = (j-1)*nsamp_win + raw.first_samp;
  to = j*nsamp_win - 1 + raw.first_samp;
  data = fiff_read_raw_segment(raw,from,to);
  E=zeros(size(data,2),nProj);
  if size(data,2) >= nsamp_win
    for k = 1:nProj
      d = data(proj_chs(k),:) - mean(data(proj_chs(k),:));%proj_chs are gradiometers
      norm_factor = norm(d);
      if norm_factor > 0
        E(:,k) = (d/norm_factor)';
      end
    end
    if exist('E','var')
      Q = orth(E);
      clear E
      P = eye(nsamp_win) - Q*Q';  % The projection operator
      data(meg_chs,:) = (P*data(meg_chs,:)')'; % Projected MEG data
    end
  end
  fiff_write_raw_buffer(outfid,data,cals);
end
fiff_finish_writing_raw(outfid);


% [xxx,times] = fiff_read_raw_segment(raw,raw.first_samp,raw.last_samp,1:306);
% %figure
% subplot(1,3,1)
% plot(times(1:3000),data(3,:),'k:')
% subplot(1,3,2)
% plot(times(1:3000),data(6,:),'k:')
% subplot(1,3,3)
% plot(times(1:3000),data(9,:),'k:')
% %fiff_finish_writing_raw(outfid);
% save('tempprojdata','data')