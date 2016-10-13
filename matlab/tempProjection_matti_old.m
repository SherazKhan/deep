function data=tempProjection_matti(data,ref,var)

%var (0,1) = amount of variation in graiometers that will be removed -
%don't put too high a value or deep brain signals will be canceled as well, and it will also lead to long computation time.
%Putting too low value will lead to low cortical signal suppression. Test a
%few values between 0-1 until you get a value that works for you.

ind_mag=3:3:306;
ind_grad=1:306;
ind_grad(ind_mag)=[];

Np=data(ind_mag,:);

Nr=data(ind_grad,:);
%Use only gradiometer closest to cortical source, chl 76
%Nr = data(76,:);


%check if there is a ref signal, if not, create one from the data
if ~exist('ref','var')
     ref = Nr;
%     %Using maximum energy channel as reference works great
%     energy = zeros(1,306);
%     for k = ind_grad
%         energy(k) = sqrt(sum(data(k,:).^2));
%     end
%     [~,indnrj] = max(energy);
%     ref = data(indnrj,:);
    
%     for i = (1:306)
%         data(i,:) = data(i,:) - mean(data(i,:));
%     end
%     [~,maxind] = max(abs(data));
%     L1gf = zeros(1,length(data));
%     for i = (1:length(data))
%         L1gf(i) = data(maxind(i),i);
%     end
%     ref = L1gf;

    data_mag_filtered = Np - Np*ref'*inv(ref*ref')*ref;
    data(ind_mag,:)=data_mag_filtered;
    
else
    if isempty(ref)
        %Shift Nr so that mean of each variable is zero.
        
%         for i = (1:length(data))
%             Nr(i) = data(maxind(i),i);
%         end
%        Nr = Nr';
        [U,S,V]=svd(Nr);%Include Np
        A=diag(S)./sum(diag(S));
        a = 0;
        i = 1;
        %Project samples (channels) onto as many PCAs so that var% of
        %variance can be explained.
        ref = [];
        while a<var%Include correlation.
            NrProj=Nr*V(:,i);
            PrincComp = sum(NrProj)*V(:,i)';
            Nr = Nr - NrProj*PrincComp;
%            NrProj = NrProj';
            a=a+A(i);
            i = i+1;
            ref = [ref; PrincComp];
%            ref = [ref; V(:,i)'];

        end
    end
    data_mag_filtered = Np - Np*ref'*inv(ref*ref')*ref;
    data(ind_mag,:)=data_mag_filtered;
    
end





