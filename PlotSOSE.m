%% Plot SOSE iteration 133 1/6 degree ws
cmp = flipud(brewermap([],'RdBu')); 
set(0,'defaultaxesFontSize',17)
set(0,'defaultAxesFontName', 'times')
set(0,'defaultTextFontName', 'times')

XC = ncread('bsose_i133_2013to2018_monthly_Wvel.nc','XC');
YC = ncread('bsose_i133_2013to2018_monthly_Wvel.nc','YC');
Zl = ncread('bsose_i133_2013to2018_monthly_Wvel.nc','Zl');
WVEL = ncread('bsose_i133_2013to2018_monthly_Wvel.nc','WVEL',[1,1,21,1],[inf,inf,1,inf]);
WVEL = nanmean(WVEL,4);

WVEL(WVEL==0) = NaN;

load SmithSandwellBathy

%%
cmp = flipud(brewermap([],'RdBu')); colormap(cmp)

h1 = figure(7); 

h11 = axes; hold on
contourf(XC,YC,WVEL',-5e-5:2e-6:5e-5,'linestyle','none')
contour(long, latss, ssbath_tenth',[-2000 -2000],'k-','linewidth',2)

ylim([-70 -35])
colormap(cmp)
caxis([-5e-5 5e-5])
cl1 = colorbar;
ylabel(cl1,'Vertical Velocity (m/s)')
xlabel('Longitude (^\circ E)')
ylabel('Latitude (^\circ N)')

box on
title('Time-averaged W at 210 m depth')

set(h1,'position',[100 100 1400 500])
set(h11,'position',[0.05 0.15 0.9 0.8],'tickdir','out')

