function data=tempProjection_matti(data,ref,Ndim,magindex,gradindex)

%var (0,1) = amount of variation in graiometers that will be removed -
%don't put too high a value or deep brain signals will be canceled as well, and it will also lead to long computation time.
%Putting too low value will lead to low cortical signal suppression. Test a
%few values between 0-1 until you get a value that works for you.


Np=data(magindex,:);

Nr=data(gradindex,:);
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
    data(magindex,:)=data_mag_filtered;
    
else
    if isempty(ref)
        
%         [Ur,Sr,Vr]=svd(Nr);
%         refr = [];
%         for i = (1:5)
%             NrProj=Nr*Vr(:,i);
%             PrincComp = sum(NrProj)*Vr(:,i)';
%             Nr = Nr - NrProj*PrincComp;
%             refr = [refr; PrincComp];
%         end
%         
%         [Up,Sp,Vp]=svd(Np);
%         refp = [];
%         for i = (1:5)
%             NpProj=Np*Vp(:,i);
%             PrincComp = sum(NpProj)*Vp(:,i)';
%             Np = Np - NpProj*PrincComp;
%             refp = [refp; PrincComp];
%         end
        
        
        
%         Nr = Nr';
%         [Ur,Sr,Vr]=svd(Nr);
%         refr = [];
%         for i = (1:5)
%             NrProj=Nr*Vr(:,i);
%             Nr = Nr - NrProj*Vr(:,i)';
%             NrProj = NrProj';
%             refr = [refr; NrProj];
%         end
%         Np = Np';
%         [Up,Sp,Vp]=svd(Np);
%         refp = [];
%         for i = (1:5)
%             NpProj=Np*Vp(:,i);
%             Np = Np - NpProj*Vp(:,i)';
%             NpProj = NpProj';
%             refp = [refp; NpProj];
%         end
%         cor = corrcoef([refp' refr']);
%         corcoffs = cor((1:5),(6:10))
%         ref = [];
%         
%         for i = (1:5)
%             for k = (1:5)
%                if corcoffs(i,k) > var
%                    ref = [ref; refr(k,:)];
%                end
%             end
%         end

        Qg = orth(Nr');
        Qm = orth(Np');
        C = Qg'*Qm;
        [Y,S,Z] = svd(C);
        theta = acos(diag(S));
        u = Qg*Y;
        v = Qm*Z;
        
        i = 1;
        ref = [];
        Nmag = Np;
        for i = (1:Ndim)
%        while theta(i) < var
            NpProj=Nmag*v(:,i);
            %signal = ifft(fft(signal)*(1-getfft(NpProj*v(:,i)'))
            Nmag = Nmag - NpProj*v(:,i)';
%            ref = [ref; sum(NpProj)*v(:,i)'];
%            i = i+1;
        end
        data_mag_filtered = Nmag;
    end
%    data_mag_filtered = Np - Np*ref'*inv(ref*ref')*ref;
    data(magindex,:)=data_mag_filtered;
    
end





