clear all;


%Set default text interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');


ind1 = load('ind1.mat');
ind1 = ind1.SNRchange;
ind2 = load('ind2.mat');
ind2 = ind2.SNRchange;
ind3 = load('ind3.mat');
ind3 = ind3.SNRchange;
ind4 = load('ind4.mat');
ind4 = ind4.SNRchange;
ind5 = load('ind5.mat');
ind5 = ind5.SNRchange;
indsurf = load('indsurf.mat');
indsurf = indsurf.SNRchange;

ind1 = sort(ind1);
ind1 = ind1(round(length(ind1)*0.1):round(length(ind1)*0.9));

ind2 = sort(ind2);
ind2 = ind2(round(length(ind2)*0.1):round(length(ind2)*0.9));

ind3 = sort(ind3);
ind3 = ind3(round(length(ind3)*0.1):round(length(ind3)*0.9));

ind4 = sort(ind4);
ind4 = ind4(round(length(ind4)*0.1):round(length(ind4)*0.9));

ind5 = sort(ind5);
ind5 = ind5(round(length(ind5)*0.1):round(length(ind5)*0.9));

indsurf = sort(indsurf);
indsurf = indsurf(round(length(indsurf)*0.1):round(length(indsurf)*0.9));

grp = [zeros(1,length(ind1)) ones(1,length(ind2)) 2*ones(1,length(ind3))...
    3*ones(1,length(ind4)) 4*ones(1,length(ind5)) 5*ones(1,length(indsurf))];
boxplot([ind1 ind2 ind3 ind4 ind5 indsurf],grp)
hold on;
xlim manual
plot([-1 10],[1 1],'r:')

ylabel(sprintf('SNR \n increase'),'rot',0)
xlabel('Dipole source')





