%clear
close all


% set axes properties
set(0,'defaultaxesFontSize',17)
set(0,'defaultAxesFontName', 'times')
set(0,'defaultTextFontName', 'times')
co = brewermap(4,'PrGn');
set(groot,'defaultAxesColorOrder',co)
cmp = flipud(brewermap([],'RdYlBu'));

% set wind stress maximum value
tau = 0.15;

% load colormap
cmapmatt = load('cmap_matt');

% set the number of iterations
countmax = 7000; %7000

% 180 years for 3dres015
% pick the right file based on wind stress value and load
% if tau == 0.15
%     %load 3dres015
%     temp = ncread('3dres015.nc','THETA');
%     ures_z = ncread('3dres015.nc','UVEL_C');
%     vres_z = ncread('3dres015.nc','VVEL_C');
%     wbt_z = ncread('3dres015.nc','WBT');
%     %wmean_z = ncread('3dres015.nc','WVEL_M');
%     wmean_z = ncread('3dres015.nc','WMEAN'); % Eulerian mean w 
%     wres_z = ncread('3dres015.nc','WVEL_C');
%     umean_z = ncread('3dres015.nc','UMEAN');
%     vmean_z = ncread('3dres015.nc','VMEAN');
%     Z = ncread('3dres015.nc','Z');
%     meanres = ncread('3dres015.nc','meanres');
    
if tau == 0.15
    load 3dres015
    temp = ncread('3dres015flat.nc','THETA');
    ures_z = ncread('3dres015flat.nc','UVEL_C');
    vres_z = ncread('3dres015flat.nc','VVEL_C');
    wbt_z = ncread('3dres015flat.nc','WBT');
    %wmean_z = ncread('3dres015flat.nc','WVEL_M');
    wmean_z = ncread('3dres015flat.nc','WMEAN');
    wres_z = ncread('3dres015flat.nc','WVEL_C');
    umean_z = ncread('3dres015flat.nc','UMEAN');
    vmean_z = ncread('3dres015flat.nc','VMEAN');
    Z = ncread('3dres015flat.nc','Z');
    meanres = ncread('3dres015flat.nc','meanres');
    
elseif tau == 0.2
    %load 3dres02
    temp = ncread('3dres02.nc','THETA');
    ures_z = ncread('3dres02.nc','UVEL_C');
    vres_z = ncread('3dres02.nc','VVEL_C');
    wbt_z = ncread('3dres02.nc','WBT');
    wmean_z = ncread('3dres02.nc','WMEAN');
    wres_z = ncread('3dres02.nc','WVEL_C');
    umean_z = ncread('3dres02.nc','UMEAN');
    vmean_z = ncread('3dres02.nc','VMEAN');
    Z = ncread('3dres02.nc','Z');
    meanres = ncread('3dres02.nc','meanres');
elseif tau == 0.1
    %load 3dflat01
    temp = ncread('3dres01.nc','THETA');
    ures_z = ncread('3dres01.nc','UVEL_C');
    vres_z = ncread('3dres01.nc','VVEL_C');
    wbt_z = ncread('3dres01.nc','WBT');
    wmean_z = ncread('3dres01.nc','WMEAN');
    wres_z = ncread('3dres01.nc','WVEL_C');
    umean_z = ncread('3dres01.nc','UMEAN');
    vmean_z = ncread('3dres01.nc','VMEAN');
    Z = ncread('3dres01.nc','Z');
    meanres = ncread('3dres01.nc','meanres');
elseif tau == 0.25
    load 3dres025
elseif tau == 0.05
    load 3dres005
else
    disp('invalid wind')
end



% analyze the overturning
meanres = meanres';
meanres = meanres(:,3:202);
resmean = mean(meanres(:,175:193),2);

% set layers grid
layers_bounds=0:8/49:8;
layer_center = layers_bounds(1:49)+diff(layers_bounds);
% get the cell split information
celldiv = layer_center(find(resmean(10:40)>0,1)+9);
uppersplit = layer_center(find(resmean==max(resmean),1));
lowersplit = layer_center(find(resmean==min(resmean(1:30)),1));


% set x-y grid
x = 10000:10000:4000000;
y = 10000:10000:2000000;

% permute the data
wres = permute(wres_z,[2 3 1]);
vres = permute(vres_z,[2 3 1]);
ures = permute(ures_z,[2 3 1]);
wmean = permute(wmean_z,[2 3 1]);
wbt = permute(wbt_z,[2 3 1]);
%wbc = permute(wbc_z,[2 3 1]);
%wgeo = permute(wgeo_z,[2 3 1]);
T = permute(temp,[2 3 1]);
wgeo = permute(EkmanMat,[2 3 1]);

% remove wall data
ures = ures(3:202,:,:);
vres = vres(3:202,:,:);
wres = wres(3:202,:,:);
wmean = wmean(3:202,:,:);
wgeo = wgeo;
 %wgeo(3:202,:,:);
wbt = wbt(3:202,:,:);
wbc = wbt; %wbc(3:202,:,:);
T = double(T(3:202,:,:));
%T(T==0) = NaN;
%T = permute(T,P);

% make grid 
[ygrid, xgrid, zgrid] = ndgrid(y,x,double(Z));




%%
% cmap = brewermap([],'RdBu');
% for i = 1:200
%     for j = 1:400
%         z25(i,j) = interp1(squeeze(T(i,j,6:25)),double(Z(6:25)),2.5);
%     end
% end
%
%
%
% figure;
% subplot(1,2,1)
% contourf(wres(:,:,5),40,'linestyle','none')
% colormap(cmap)
% caxis([-5e-4 5e-4])
% colorbar
%
%
% subplot(1,2,2)
% contourf(wmean(:,:,5),40,'linestyle','none')
% caxis([-5e-4 5e-4])
% colorbar
%% notes to self
% renormalize so that we multiply by the mean time for all trajectories so
% that displacement is the right units

% to create means for each cell, create way to weight by volume of each
% temperature layer, start by calculating average of each volume layer

% must use the temperature grid to define what the mean temperatures are

% maybe work on making these statistics more rigorous

%%
disp('started')

% 4 is upwelling on the upper branch
% 2 is upwelling on the lower branch
% 1 us sinking on the lower branch
% 6 is sinking on the upper branch

% make layer_center the temperature that we release particles on
temps = layer_center;

% make z grid of the temperature locations
zmap = nan(200,length(temps));
for i = 2:200
    zmap(i,:) = interp1(squeeze(T(i,6:26,3)),double(Z(6:26)),temps);
end

% create x,y,z locations to release the particles
xfloats = [];
yfloats = [];
zfloats = [];
for k = 1:length(temps)
    xfloats = [xfloats;xgrid(1:end,3,1)];
    yfloats = [yfloats;ygrid(1:end,3,1)];
    zfloats = [zfloats;zmap(1:end,k)];
end



%     h32 = figure(43);
%     s = surf(xgrid(:,:,1)/1000,ygrid(:,:,1)/1000,zmap,'FaceAlpha',0.5);
%     s.EdgeColor = 'none';
%     hold on




% run particle tracking
[posx,posy,posz,zmean,zbc,zbt,zgeo, meanx, resdispz, meandispz, bcdispz, btdispz, geodispz,timemean] = ParticleTracking(xfloats,yfloats,...
    zfloats,ures,vres,wres,wmean,wgeo,wbc,wbt,xgrid,ygrid,zgrid,countmax);


%%

% for k = 1:length(temps)
%     hh2 = figure(33);
%     ti = sprintf('Residual Displacement Map %g C tau = %g Nm^2',temps(k),tau);
%     dispz = posz - posz(:,1);
%     scatter(posx(:)/1000,posy(:)/1000,3,dispz(:))
%     xlim([0 4000])
%     ylim([0 2000])
%     title(ti)
%     c1 = colorbar;
%     colormap(cmapmatt.cmap)
%     xlabel('X (km)')
%     ylabel('Y (km)')
%     ylabel(c1,'Displacement (m)')
%     grid on
%     caxis([-400 400])
%     SaveFigureThesis(hh2,ti)
%     
%     close
%     
% %     hh4 = figure(34);
% %     ti = sprintf('Mean Displacement Map %g C tau = %g Nm^2',temps(k),tau);
% %     %dispzmean = poszmean - poszmean(:,1);
% %     scatter(posx(:)/1000,posy(:)/1000,3,zmean(:))
% %     xlim([0 4000])
% %     ylim([0 2000])
% %     title(ti)
% %     c1 = colorbar;
% %     colormap(cmapmatt.cmap)
% %     xlabel('X (km)')
% %     ylabel('Y (km)')
% %     ylabel(c1,'Displacement (m)')
% %     caxis([-400 400])
% %     grid on
% %     SaveFigureThesis(hh4,ti)
% %     
% %     close
% %     hh4 = figure(34);
% %     ti = sprintf('Advection-Mean Displacement Map %g C tau = %g Nm^2',temps(k),tau);
% %     %dispzmean = poszmean - poszmean(:,1);
% %     scatter(posx(:)/1000,posy(:)/1000,3,zamean(:))
% %     xlim([0 4000])
% %     ylim([0 2000])
% %     title(ti)
% %     c1 = colorbar;
% %     colormap(cmapmatt.cmap)
% %     xlabel('X (km)')
% %     ylabel('Y (km)')
% %     ylabel(c1,'Displacement (m)')
% %     caxis([-400 400])
% %     grid on
% %     SaveFigureThesis(hh4,ti)
% %     
% %     close
% %     hh4 = figure(34);
% %     ti = sprintf('Thickness-Eddy Displacement Map %g C tau = %g Nm^2',temps(k),tau);
% %     %dispzmean = poszmean - poszmean(:,1);
% %     scatter(posx(:)/1000,posy(:)/1000,3,zeddy(:))
% %     xlim([0 4000])
% %     ylim([0 2000])
% %     title(ti)
% %     c1 = colorbar;
% %     colormap(cmapmatt.cmap)
% %     xlabel('X (km)')
% %     ylabel('Y (km)')
% %     ylabel(c1,'Displacement (m)')
% %     caxis([-400 400])
% %     grid on
% %     SaveFigureThesis(hh4,ti)
% %     
% %     close
% %     
% %     hhh1 = figure(23);
% %     hold on
% %     plot(meanx/1000,resdispz,'-')
% %     plot(meanx/1000,meandispz,'--')
% %     plot(meanx/1000,meanadispz,':')
% %     plot(meanx/1000,eddydispz,'-.')
% %     
% %     meanxtot = meanx/1000;
% %     
% %     
% %     legend('Residual','Mean','Mean-Advection','Thickness')
% %     grid on
% end

%%

% legend('1 C','1 C - Mean','2 C','2 C - Mean','3.5 C','3.5 C - Mean','5.3 C','5.3 C - Mean')
% xlabel('X (km)')
% ylabel('Dispacement (m)')
% ti = sprintf('Displacement tau = %g Nm^2',tau);
% title(ti)
% SaveFigureThesis(hhh1,ti)
% 
% %%
% close all
% hhh1 = figure(23);
% hold on
% for k = 1:2
%     plot(meanx/1000,resdispz(:,k),'-','color',co(k,:),'linewidth',2)
%     plot(meanx/1000,meandispz(:,k),'--','color',co(k,:),'linewidth',2)
%     plot(meanx/1000,eddydispz(:,k),':','color',co(k,:),'linewidth',2)
%     plot(meanx/1000,meanadispz(:,k),'-.','color',co(k,:),'linewidth',2)
% end
% 
% legend('Downwelling','Downwelling - Mean','Downwelling - Eddy',...
%     'Downwelling - Advection','Upwelling','Upwelling - Mean',...
%     'Upwelling - Eddy','Upwelling - Advection')
% xlabel('X (km)')
% ylabel('Dispacement (m)')
% ti = sprintf('Lower Cell Displacement tau = %g Nm^2',tau);
% title(ti)
% grid
% ylim([-150 350])
% 
% SaveFigureThesis(hhh1,ti)
% 
% %%
% close all
% hhh1 = figure(23);
% hold on
% for k = 4:-1:3
%     plot(meanx/1000,resdispz(:,k),'-','color',co(k,:),'linewidth',2)
%     plot(meanx/1000,meandispz(:,k),'--','color',co(k,:),'linewidth',2)
%     plot(meanx/1000,eddydispz(:,k),':','color',co(k,:),'linewidth',2)
%     plot(meanx/1000,meanadispz(:,k),'-.','color',co(k,:),'linewidth',2)
% end
% 
% legend('Downwelling','Downwelling - Mean','Downwelling - Eddy',...
%     'Downwelling - Advection','Upwelling','Upwelling - Mean',...
%     'Upwelling - Eddy','Upwelling - Advection')
% xlabel('X (km)')
% ylabel('Dispacement (m)')
% ti = sprintf('Upper Cell Displacement tau = %g Nm^2',tau);
% title(ti)
% grid
% xlim([0 4000])
% %ylim([-150 350])
% SaveFigureThesis(hhh1,ti)

%%
h11 = figure(88);

contourf(meanx/1000,temps,resdispz',50,'linestyle','none')
hold on
plot(meanx/1000,ones(1,length(meanx))*celldiv,'-k','linewidth',2)
plot(meanx/1000,ones(1,length(meanx))*uppersplit,'--k','linewidth',2)
plot(meanx/1000,ones(1,length(meanx))*lowersplit,'--k','linewidth',2)
caxis([-2.5e-6 2.5e-6])
ylim([0 5])
xlabel('X (km)')
ylabel('Temperature Layers (C)')
ti = sprintf('Residual Displacement tau = %g Nm^2',tau);
title(ti)
c1 = colorbar;
colormap(cmapmatt.cmap)
ylabel(c1,'Displacement/time (m/s)')
SaveFigureThesis(h11,ti)







h11 = figure(89);
colormap(cmapmatt.cmap)
contourf(meanx/1000,temps,meandispz',50,'linestyle','none')
hold on
plot(meanx/1000,ones(1,length(meanx))*celldiv,'-k','linewidth',2)
plot(meanx/1000,ones(1,length(meanx))*uppersplit,'--k','linewidth',2)
plot(meanx/1000,ones(1,length(meanx))*lowersplit,'--k','linewidth',2)
caxis([-2.5e-6 2.5e-6])
xlabel('X (km)')
ylim([0 5])
ylabel('Temperature Layers (C)')
ti = sprintf('Mean Displacement tau = %g Nm^2',tau);
title(ti)
c1 = colorbar;
colormap(cmapmatt.cmap)
ylabel(c1,'Displacement/time (m/s)')


SaveFigureThesis(h11,ti)

h11 = figure(90);
colormap(cmapmatt.cmap)
contourf(meanx/1000,temps,resdispz'-meandispz',50,'linestyle','none')
hold on
plot(meanx/1000,ones(1,length(meanx))*celldiv,'-k','linewidth',2)
plot(meanx/1000,ones(1,length(meanx))*uppersplit,'--k','linewidth',2)
plot(meanx/1000,ones(1,length(meanx))*lowersplit,'--k','linewidth',2)
caxis([-2.5e-6 2.5e-6])
xlabel('X (km)')
ylabel('Temperature Layers (C)')
ti = sprintf('Eddy Displacement tau = %g Nm^2',tau);
title(ti)
c1 = colorbar;
colormap(cmapmatt.cmap)
ylabel(c1,'Displacement/time (m/s)')
ylim([0 5])


SaveFigureThesis(h11,ti)

% h11 = figure(91);
% colormap(cmapmatt.cmap)
% contourf(meanx/1000,temps,bcdispz'+btdispz',50,'linestyle','none')
% hold on
% plot(meanx/1000,ones(1,length(meanx))*celldiv,'-k','linewidth',2)
% plot(meanx/1000,ones(1,length(meanx))*uppersplit,'--k','linewidth',2)
% plot(meanx/1000,ones(1,length(meanx))*lowersplit,'--k','linewidth',2)
% xlabel('X (km)')
% ylabel('Temperature Layers (C)')
% ti = sprintf('Eddy Displacement tau = %g Nm^2',tau);
% title(ti)
% c1 = colorbar;
% caxis([-2.5e-6 2.5e-6])
% colormap(cmapmatt.cmap)
% ylabel(c1,'Displacement/time (m/s)')
% 
% 
% SaveFigureThesis(h11,ti)

% 
% 
h11 = figure(92);
colormap(cmapmatt.cmap)
contourf(meanx/1000,temps,geodispz',50,'linestyle','none')
hold on
plot(meanx/1000,ones(1,length(meanx))*celldiv,'-k','linewidth',2)
plot(meanx/1000,ones(1,length(meanx))*uppersplit,'--k','linewidth',2)
plot(meanx/1000,ones(1,length(meanx))*lowersplit,'--k','linewidth',2)
xlabel('X (km)')
ylabel('Temperature Layers (C)')
ti = sprintf('Ekman Displacement tau = %g Nm^2',tau);
title(ti)
c1 = colorbar;
caxis([-2.5e-6 2.5e-6])
ylim([0 5])
colormap(cmapmatt.cmap)
ylabel(c1,'Displacement/time (m/s)')
SaveFigureThesis(h11,ti)
% 
% 
h11 = figure(93);
colormap(cmapmatt.cmap)
contourf(meanx/1000,temps,meandispz'-geodispz',50,'linestyle','none')
hold on
plot(meanx/1000,ones(1,length(meanx))*celldiv,'-k','linewidth',2)
plot(meanx/1000,ones(1,length(meanx))*uppersplit,'--k','linewidth',2)
plot(meanx/1000,ones(1,length(meanx))*lowersplit,'--k','linewidth',2)
xlabel('X (km)')
ylabel('Temperature Layers (C)')
ti = sprintf('Non-Ekman Mean Displacement tau = %g Nm^2',tau);
title(ti)
c1 = colorbar;
caxis([-2.5e-6 2.5e-6])
ylim([0 5])
colormap(cmapmatt.cmap)
SaveFigureThesis(h11,ti)
ylabel(c1,'Displacement/time (m/s)')

%% average the values in each layer
uppercellinds = find(temps>celldiv);
lowerup = find(temps<celldiv & temps>lowersplit);
lowerdown = find(temps<lowersplit);

h11 = figure(94);
hold on
plot(meanx/1000,nanmean(resdispz(:,uppercellinds),2),'-g','linewidth',2)
plot(meanx/1000,nanmean(geodispz(:,uppercellinds),2),'--k','linewidth',2)
plot(meanx/1000,nanmean(bcdispz(:,uppercellinds),2),'-.b','linewidth',2)
plot(meanx/1000,nanmean(btdispz(:,uppercellinds),2),'-.m','linewidth',3)
plot(meanx/1000,nanmean(meandispz(:,uppercellinds),2),':c','linewidth',2)
legend('Residual','Geostrophic','Eddy Buoyancy Fluxes','Mean Advection','Eddy Momentum Fluxes','location','best')

hold on
ylim([-2e-6 2e-6])
grid on
xlabel('X (km)')
%ylabel('Temperature Layers (C)')
ti = sprintf('Upper Cell Displacement tau = %g Nm^2',tau);
title(ti)
%c1 = colorbar;
%colormap(cmapmatt.cmap)
ylabel('Displacement/time (m/s)')
SaveFigureThesis(h11,ti)


h11 = figure(95);
hold on
plot(meanx/1000,nanmean(resdispz(:,lowerup),2),'-g','linewidth',2)
plot(meanx/1000,nanmean(geodispz(:,lowerup),2),'--k','linewidth',2)
plot(meanx/1000,nanmean(eddydispz(:,lowerup),2),'-.b','linewidth',2)
plot(meanx/1000,nanmean(meanadispz(:,lowerup),2),'-.m','linewidth',3)
plot(meanx/1000,nanmean(resdispz(:,lowerup),2)-nanmean(eddydispz(:,lowerup),2)-nanmean(geodispz(:,lowerup),2)-nanmean(meanadispz(:,lowerup),2),':c','linewidth',2)
legend('Residual','Geostrophic','Baroclinic Eddies','Mean Advection','Barotropic Eddies','location','best')

hold on
ylim([-2e-6 2e-6])
grid on
xlabel('X (km)')
%ylabel('Temperature Layers (C)')
ti = sprintf('Lower Cell Upwelling Displacement tau = %g Nm^2',tau);
title(ti)
%c1 = colorbar;
%colormap(cmapmatt.cmap)
ylabel('Displacement/time (m/s)')
SaveFigureThesis(h11,ti)


h11 = figure(96);
hold on
plot(meanx/1000,nanmean(resdispz(:,lowerdown),2),'-k')
plot(meanx/1000,nanmean(geodispz(:,lowerdown),2),'--k')
plot(meanx/1000,nanmean(eddydispz(:,lowerdown),2),'-.k')
plot(meanx/1000,nanmean(meanadispz(:,lowerdown),2),'-.k','linewidth',2)
plot(meanx/1000,nanmean(resdispz(:,lowerdown),2)-nanmean(eddydispz(:,lowerdown),2)-nanmean(geodispz(:,lowerdown),2)-nanmean(meanadispz(:,lowerdown),2),':k')
legend('Residual','Geostrophic','Baroclinic Eddies','Mean Advection','Barotropic Eddies')

hold on
ylim([-2e-6 2e-6])
grid on
xlabel('X (km)')
%ylabel('Temperature Layers (C)')
ti = sprintf('Lower Cell Downwelling Displacement tau = %g Nm^2',tau);
title(ti)
%c1 = colorbar;
%colormap(cmapmatt.cmap)
ylabel('Displacement/time (m/s)')
SaveFigureThesis(h11,ti)