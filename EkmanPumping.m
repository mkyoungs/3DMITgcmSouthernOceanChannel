%% plot the new figure 9
%% headers and read in code

clear
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



% 180 years for 3dres015
% pick the right file based on wind stress value and load
if tau == 0.15
    %load 3dres015
    temp = ncread('3dres015.nc','THETA');
    Z = ncread('3dres015.nc','Z');
    wmean = ncread('3dres015.nc','WMEAN');
else
    disp('None')
end


%% calculated intercepts for northern boundary and surface using contour

surf2 = contourc(temp(:,3:202,1),[2 2]);
surf2 = surf2(:,2:end);
%surf1 = contourc(temp(:,3:202,1),[1.5 1.5]);
surf4 = contourc(temp(:,3:202,1),[4 4]);
surf4 = surf4(:,2:end);


% allsurfaces = contourc(temp(:,3:202,1),0.5:0.2:4.6);
% allsurfaces = allsurfaces(:,2:end);
% contour(allsurfaces)

%nbdy1 = contourc(squeeze(temp(:,202,:)),[1.5 1.5]);
nbdy2 = contourc(Z,1:400,squeeze(temp(:,202,:)),[2 2]);
nbdy4 = contourc(Z,1:400,squeeze(temp(:,202,:)),[4 4]);
nbdy2 = nbdy2(:,2:end);
nbdy4 = nbdy4(:,2:end);

count = 0;
while max(diff(surf2(2,:)))> 0 || max(diff(surf4(2,:)))> 0
    count = count+1;
    ind2 = find(diff(surf2(2,:))<0);
    ind4 = find(diff(surf4(2,:))<0);
    surf2 = surf2(:,ind2);
    surf4 = surf4(:,ind4);

end

% interpolate onto xgrid
surfNew2 = interp1(surf2(2,:),surf2(1,:),1:400);
surfNew4 = interp1(surf4(2,:),surf4(1,:),1:400);

nbdyNew2 =  interp1(nbdy2(2,:),nbdy2(1,:),1:400);
nbdyNew4 = interp1(nbdy4(2,:),nbdy4(1,:),1:400);


%% calculate and plot slope
slope2 = nbdyNew2./(2000000-surfNew2*10000);
slope4 = nbdyNew4./(2000000-surfNew4*10000);

h1 = figure(33); 
hh1 = axes;
plot(10:10:4000,slope2,'k-',10:10:4000,slope4,'k:','linewidth',2)
legend('2 C','4 C','location','best')
grid on
box on
text(100,-2e-4,'(a)','fontsize',16)
ylabel('Slope')

%% calculate ekman pumping
nytot = 204;
wind = sin((0:nytot-1)*pi/(nytot-1))*0.15;
dwinddy = pi/2030000*cos((0:nytot-1)*pi/(nytot-1))*0.15;
beta = 1e-11;
ftot = -1e-4 + beta*(1:204)*10000;
rho0 = 1000;
ekmanpumping = -1/rho0./ftot./ftot.*(dwinddy.*ftot - wind*beta);
ekmanpumpingdy = diffxy(1:204,ekmanpumping);


%winddecay(i,:) = [sin((0:nytot-1)*pi/(nytot-1))*0.1].*(tanh((0.75-ygrid(:,1)')/.05)+1)/2;



%% calculate w temp using above pieces
ekman2 = interp1(1:200,ekmanpumping(3:end-2),surfNew2);
ekman4 = interp1(1:200,ekmanpumping(3:end-2),surfNew4);


%% plot the ekman pumping

hh2 = axes;
plot(10:10:4000,ekman2,'k-',10:10:4000,ekman4,'k:','linewidth',2);
xlabel('X (km)')
ylabel('Ekman Pumping (m s^{-1})')
grid on
box on
text(100,2.5e-6,'(b)','fontsize',16)



%% shape the figure
set(h1,'position',[100 100 700 600])
set(hh1, 'position',[0.1 0.56 0.8 0.38])
set(hh2, 'position',[0.1 0.1 0.8 0.38])

%%
for i = 1:400
    i
    ytemps = interp1(smoothdata(squeeze(temp(i,3:202,1)),'movmean',3),1:200,0.1:0.1:8);
    Wek = interp1(1:200,ekmanpumping(3:202),ytemps);
    %WekdT = diffxy(0.4:0.1:4.9,Wek);
    wvec = interp1(0.1:0.1:8,Wek,reshape(temp(i,3:202,:),1,[]));
    EkmanMat(i,:,:) = reshape(wvec,200,32);
end

%%
EkmanMat = ekmanpumping(3:202)+zeros(400,200,32);