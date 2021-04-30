function [ output_args ] = MixedLayer( HF, theta, psiEm)
%MIXEDLAYER Summary of this function goes here
%   Detailed explanation goes here
global yg
l = size(HF,3);
HF = squeeze(nanmean(HF,2));
theta = squeeze(nanmean(theta,2));

for i = 1:l
theta(:,i) = smooth(theta(:,i),120);
HF(:,i) = smooth(HF(:,i),120);
end
psiEm = squeeze(nanmean(psiEm,1));

%dpsiEk = 0.1*sin((1:204)*pi/204)*10;
dbdy = diffxy(yg(1:l,3:202)',theta(3:202,:));

for i = 2:l
dHF(:,i) = (HF(:,i)-HF(:,1))/10^3/4000;
dbdt(:,i) = (theta(:,i)-theta(:,i-1))/2.63e6;
dpsiEm(:,i) = psiEm(:,i)-psiEm(:,1);
end
dpsiML = (dHF(3:202,:) - 100*dbdt(3:202,:))./dbdy;
%dpsiML(:,2) = dpsiEm(3:202,2);
dpsiEddy = dpsiML-dpsiEm(3:202,:);

h43 = figure;
contourf((1:l)/12,1:200,dpsiEddy)
colorbar
xlabel('Year')
ylabel('Y')
title('d\psi Eddy')
SaveFigureThesis(h43,'dpsi Eddy')

h44 = figure;
contourf((1:l)/12,1:204,dHF)
colorbar
xlabel('Year')
ylabel('Y')
title('dHeat Flux')
SaveFigureThesis(h44,'d HeatFlux')

h45 = figure;
contourf((1:l)/12,1:204,dbdt)
colorbar
xlabel('Year')
ylabel('Y')
title('dbdt')
SaveFigureThesis(h45,'dbdt')

h46 = figure;
contourf((1:l)/12,1:200,dpsiML)
colorbar
xlabel('Year')
ylabel('Y')
title('dpsi ML')
SaveFigureThesis(h46,'d psi Mixed Layer')



end

