function [posx,posy,posz,poszmean,poszbc,poszbt,poszg,meanx, resdispz, meandispz, bcdispz, btdispz, geodispz,timemean] = ParticleTracking(xfloats,yfloats,zfloats,u,v,w,wmean,wgeo,wbc,wbt,xgrid,ygrid,zgrid,countmax,vcodec)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% set time interval
dt = 100000;

% initialize count and time
count = 0;
t = 0;

% permute data again
P = [1 3 2];
% xgrid = permute(xgrid, P);
% ygrid = permute(ygrid, P);
% zgrid = permute(zgrid, P);
u = permute(u, P);
v = permute(v, P);
w = permute(w, P);
wmean = permute(wmean, P);
wbc = permute(wbc,P);
wbt = permute(wbt,P);
wgeo = permute(wgeo,P);

% P = [3 1 2];
% xgrid = permute(xgrid, P);
% ygrid = permute(ygrid, P);
% zgrid = permute(zgrid, P);

% initialize x,y,z positions as empty
posx = NaN(length(xfloats),countmax);
posy = posx;
posz = posx;
poszmean = posx;
poszbt = posx;
poszbc = posx;
poszgeo = posx;
zmean = zeros(length(xfloats),1);
zbt = zmean;
zbc= zmean;
zgeo = zmean;



% put in the initial positions
posx(:,1) = xfloats;
posy(:,1) = yfloats;
posz(:,1) = zfloats;
poszmean(:,1) = zmean;
poszbt(:,1) = zbt;
poszbc(:,1) = zbc;
poszg(:,1) = zmean;


% start the iteration
while count < countmax
    % add 1 to the count
    count = count+1;
    
    % only display count if it is a factor of 10
    if rem(count,10)==0
        disp(count)
    end
    
    % calculate velocity at the float positions
    ufloats=interpn(ygrid,xgrid,zgrid,u,yfloats,xfloats,zfloats);
    vfloats=interpn(ygrid,xgrid,zgrid,v,yfloats,xfloats,zfloats);
    wfloats=interpn(ygrid,xgrid,zgrid,w,yfloats,xfloats,zfloats);
    
    %     ufloats2=interp3(xgrid,ygrid,zgrid,u,...
    %         xfloats+dt*ufloats,yfloats+dt*vfloats,zfloats+dt*wfloats);
    %     vfloats2=interp3(xgrid,ygrid,zgrid,v,...
    %         xfloats+dt*ufloats,yfloats+dt*vfloats,zfloats+dt*wfloats);
    %     wfloats2=interp3(xgrid,ygrid,zgrid,w,...
    %         xfloats+dt*ufloats,yfloats+dt*vfloats,zfloats+dt*wfloats);
    %     xfloats=xfloats+(dt/2)*(ufloats+ufloats2);
    %     yfloats=yfloats+(dt/2)*(vfloats+vfloats2);
    %     zfloats=zfloats+(dt/2)*(wfloats+wfloats2);
    
    % step the floats forward in time
    xfloats=xfloats+(dt)*(ufloats);
    yfloats=yfloats+(dt)*(vfloats);
    zfloats=zfloats+(dt)*(wfloats);
    
    % write in mean and eddy components
    wfloatmean=interpn(ygrid,xgrid,zgrid,wmean,yfloats,xfloats,zfloats);
    wfloatbt=interpn(ygrid,xgrid,zgrid,wbt,yfloats,xfloats,zfloats);
    wfloatbc=interpn(ygrid,xgrid,zgrid,wbc,yfloats,xfloats,zfloats);
    wfloatgeo=interpn(ygrid,xgrid,zgrid,wgeo,yfloats,xfloats,zfloats);
    zmean = zmean + dt*wfloatmean;
    zbc = zbc + dt*wfloatbc;
    zbt = zbt + dt*wfloatbt;
    zgeo = zgeo + dt*wfloatgeo;
    
    
    
    
    
    %     if sum(xfloats >= 3990000)>=1
    %             posx(xfloats>=3990000,count-1)= NaN;
    %             xfloats(xfloats>=3990000) = xfloats(xfloats>=3990000)-3980000;
    %     end
    
    % find indices where the the floats are out of bounds or not moving
    kind = find(abs(ufloats) < 1e-3 & abs(vfloats) < 1e-3 & abs(wfloats) < 1e-3);
    
    jind = find(abs(ufloats) < 1e-3 & ufloats<0 & abs(vfloats) < 1e-3 & abs(wfloats) < 1e-3);
    
    lind = find(yfloats>1930000);
    if ~isempty(kind)||~isempty(jind)
        %disp('delete')
    end
    
    % delete the out of bounds locations
    posx(jind,:) = NaN;
    posx(kind,:) = NaN;
    zmean(kind,:) = NaN;
    posx(jind,:) = NaN;
    posx(lind,:) = NaN;
    zmean(lind,:) = NaN;
    zbc(kind,:) = NaN;
    zbc(lind,:) = NaN;
    zmean(jind,:) = NaN;
    zbc(jind,:) = NaN;
    zbt(kind,:) = NaN;
    zbt(lind,:) = NaN;
    zbt(jind,:) = NaN;
    xfloats(kind) = NaN;
    xfloats(jind) = NaN;
    yfloats(kind) = NaN;
    zfloats(kind) = NaN;
    xfloats(lind) = NaN;
    zgeo(jind,:) = NaN;
    zgeo(kind,:) = NaN;
    zgeo(lind,:) = NaN;

    % delete floats if they enter mixed layer
    if sum(zfloats>-250)>=1
        posz(zfloats>-250,:) = NaN;
        poszmean(zfloats>-250,:) = NaN;
        poszbc(zfloats>-250,:) = NaN;
        poszbt(zfloats>-250,:) = NaN;
        poszg(zfloats>-250,:) = NaN;
        zfloats(zfloats>-250) = NaN;
    end
    
    % add the positions to output
    posx(:,count+1) = xfloats;
    posy(:,count+1) = yfloats;
    posz(:,count+1) = zfloats;
    poszmean(:,count+1) = zmean;
    poszbt(:,count+1) = zbt;
    poszbc(:,count+1) = zbc;
    poszg(:,count+1) = zgeo;
    
    % make movie
    if exist('vcodec','var') == 1
        plot3(xfloats/1000,yfloats/1000,zfloats,'k.')
        view(15,56);
        
        xlabel('X (km)')
        ylabel('Y (km)')
        zlabel('Z (m)')
        drawnow
        ax = gca;
        ax.Units = 'pixels';
        pos = ax.Position;
        ti = ax.TightInset;
        rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
        
        
        frame = getframe(ax,rect);
        writeVideo(vcodec,frame);
        
    end
    
    % if all the floats are dead, stop iterating
    if sum(isnan(xfloats)) == length(xfloats)
        break
    end
    
    % add the time
    t=t+dt;
end


%for j = 1:length(xfloats)
%end

% initialize x positions
meanx = 0:10000:4000000;

% calculate displacement
dispz = posz - posz(:,1);

% initialize the interpolation
zdisp = nan(length(posz(:,1)),length(meanx));
zmeandisp = zdisp;
zbcdisp = zdisp;
zbtdisp = zdisp;
zgeodisp = zdisp;
time = ones(1,length(xfloats));

% iterate over all the floats
for i = 1:length(posz(:,1))
    % eliminate the end wehre the floats don't move
    ind0 = find(diff(posx(i,:))==0,1);
    posx(i,ind0:end) = NaN;
    indnan = find(isnan(posx(i,:)),1);
    if indnan < 3
        
    elseif ~isempty(indnan)
        % interpolate each float on regular x grid and divide by total time
        time(i) = dt*find(isnan(posx(i,:)),1,'first');
        zdisp(i,:) = interp1(posx(i,1:indnan-1),dispz(i,1:indnan-1),meanx)/time(i);
        zmeandisp(i,:) = interp1(posx(i,1:indnan-1),poszmean(i,1:indnan-1),meanx)/time(i);
        zbtdisp(i,:) = interp1(posx(i,1:indnan-1),poszbt(i,1:indnan-1),meanx)/time(i);
        zbcdisp(i,:) = interp1(posx(i,1:indnan-1),poszbc(i,1:indnan-1),meanx)/time(i);
        zgeodisp(i,:) = interp1(posx(i,1:indnan-1),poszg(i,1:indnan-1),meanx)/time(i);
    else 
        % interpolate float on regular x grid and divide by final time
        time(i) = t;
        zdisp(i,:) = interp1(posx(i,:),dispz(i,:),meanx)/t;
        zmeandisp(i,:) = interp1(posx(i,:),poszmean(i,:),meanx)/t;
        zbtdisp(i,:) = interp1(posx(i,:),poszbt(i,:),meanx)/t;
        zbcdisp(i,:) = interp1(posx(i,:),poszbc(i,:),meanx)/t;
        zgeodisp(i,:) = interp1(posx(i,:),poszg(i,:),meanx)/t;
    end
end


for k = 1:length(xfloats)/200
    % average over all the floats
    resdispz(:,k) =  nanmean(zdisp(1+(k-1)*200:k*200,:),1);
    meandispz(:,k) = nanmean(zmeandisp(1+(k-1)*200:k*200,:),1);
    btdispz(:,k) = nanmean(zbtdisp(1+(k-1)*200:k*200,:),1);
    bcdispz(:,k) = nanmean(zbcdisp(1+(k-1)*200:k*200,:),1);
    geodispz(:,k) = nanmean(zgeodisp(1+(k-1)*200:k*200,:),1);
   
end
% calculate the mean time
timemean = mean(time);


end

