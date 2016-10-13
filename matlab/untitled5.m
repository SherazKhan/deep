Np = [1 1 0; 1 0 1];
Nr = [0 1 0];

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
for i = (1:1)
    %        while theta(i) < var
    NpProj=Nmag*v(:,i);
    Nmag = Nmag - NpProj*v(:,i)';
    %            ref = [ref; sum(NpProj)*v(:,i)'];
    %            i = i+1;
end

data_mag_filtered = Nmag;
