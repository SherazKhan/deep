function data_n=temporalProjection_matrix(data,window_length,Fs)


proj_chs = 74;
mag_chs = 3:3:306;



nsamp = size(data,2);
data_length = nsamp/Fs;




nwin = round(data_length/window_length);  % Number of time windows
% for the temporal projection
nsamp_win = round(nsamp/nwin);  % Number of samples in one window



nProj = length(proj_chs);

data_n=zeros(size(data));

for j = 1:nwin
    from = floor((j-1)*nsamp_win + 1);
    to = ceil(j*nsamp_win);
    tmp = data(:,from:to);
    E=zeros(size(tmp,2),nProj);
    if size(tmp,2) >= nsamp_win
        for k = 1:nProj
            d = tmp(proj_chs(k),:) - mean(tmp(proj_chs(k),:));%proj_chs are gradiometers
            norm_factor = norm(d);
            if norm_factor > 0
                E(:,k) = (d/norm_factor)';
            end
        end
        if exist('E','var')
            Q = orth(E);
            P = eye(nsamp_win) - Q*Q';  % The projection operator
            data_n(mag_chs,from:to) = (P*tmp(mag_chs,:)')'; % Projected MEG data
        end
    end
    
end

proj_chs = 1:306;
proj_chs(3:3:306)=[];
data_n(proj_chs,:) = data(proj_chs,:);

data_n=data_n(:,1:to);

