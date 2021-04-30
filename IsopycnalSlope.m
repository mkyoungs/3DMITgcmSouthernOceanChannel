function [slopeparam, growthrateparam] = IsopycnalSlope(theta,U,V)
% Script to compute isopycnal slope

global xc yc yg zc zf drc dyc drf ny nx

%input theta, dyc, drc
thetazonal = squeeze(nanmean(theta(:,:,:),1));

%[~,DYCgrid,~] = meshgrid(xc,dyc(2:end-1),zc);
dbdy = [thetazonal(2:end,:)-thetazonal(1:end-1,:)]./dyc(1,1);

dbdyng = interp2(squeeze(zc),yg(1,2:end),dbdy,squeeze(zc),yc(1,:));


[DRCgrid,~] = meshgrid(drc(2:end-1),xc(1,:));
dbdz = [thetazonal(:,2:end)-thetazonal(:,1:end-1)]./DRCgrid;

dbdzng = interp2(squeeze(zc(2:end)+drc(2:end-1)/2),yc(1,:),dbdz,squeeze(zc),yc(1,:));

slope = dbdyng./dbdzng;
%slope(slope > 0) = NaN;

slopeparam = nansum(nanmean(slope(3:ny-2,3:end),1).*squeeze(drf(3:end))')/3750;
slopeN = nansum(nanmean(slope(ny/2+1:ny-2,3:end),1).*squeeze(drf(3:end))')/3750;
slopeS =  nansum(nanmean(slope(3:ny/2,3:end),1).*squeeze(drf(3:end))')/3750;




%% compute N^2 on z locations


f = meshgrid(-1e-4 + 1e-11*yc(1,:),xc(:,1),1:32);
alpha = 2e-4;
thetaprofile = nanmean(reshape(theta(:,3:ny-2,:),nx*(ny-4),32),1);
[~,~,thetamat] = meshgrid(ones(ny,1),ones(nx,1),thetaprofile);



%N2p = diffxy(zc,9.8*alpha*thetaprofile);



N2 = diffxy(zc,9.8*alpha*thetamat,3);

%N2 = N2mat;



%% compute eady growth rate after Williams et al. 2007
%layer = 1;
omega = 0.31*abs(f./sqrt(N2).*sqrt(diffxy(zc,U,3).^2 +diffxy(zc,V,3).^2));

%T = omega*3600*24;

%rate = T(:,:,layer);
growthrate = -trapz(squeeze(zc),omega,3)/4000;
%rate = trapz(level(1:26),T(:,:,1:26),3)/level(26);
%maxrate = max(rate(:));

% h1 = figure; hold on;
% contourf(rate','linestyle','none')
% 
% %caxis([0 maxrate])
% colorbar
% colormap(brewermap([],'BuPu'))
% ti = 'Vertically Averaged Eady Growth Rate (days^{-1})';
% title(ti)
% contour(theta(:,:,layer)',1:20,'k')
% 
% 
% growthrate = rate(450:550,60:180);
 growthrateparam = mean(growthrate(:));


%disp(slopeparam)
end
