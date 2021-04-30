%% compute overturning

cd ~/'Dropbox (MIT)'/Work/Thesis/MITgcm/test
theta1 = rdmds('THETA',NaN);
theta1 = nanmean(theta1,4);

u = rdmds('UVEL',NaN);
u = nanmean(u,4);

v = rdmds('VVELMASS',NaN);
v = nanmean(v,4);

u2 = rdmds('UVELSQ',NaN);
u2 = nanmean(u2,4);

v2 = rdmds('VVELSQ',NaN);
v2 = nanmean(v2,4);

xc=rdmds('XC');
xg=rdmds('XG');

dxc=rdmds('DXC');
dxg=rdmds('DXG');

yc=rdmds('YC');
yg=rdmds('YG');

dyc=rdmds('DYC');
dyg=rdmds('DYG');

rac = rdmds('RAC');

zc=rdmds('RC');
zf=rdmds('RF');
drf = rdmds('DRF');

topo = rdmds('Depth');

[stat,mess] = fileattrib('T.0*');
nfiles = size(mess,2);
nz = length(zc);
[nx, ny] = size(xc);

thetay= zeros(ny,nz+1);
psiEM= zeros(ny,nz+1);
Psi_EM = zeros(nx,ny,nz+1);

Psi_EM(:,:,:) = 0;
psi_EM(:,:,nz+1)=0;
for kk=nz:-1:1
    Psi_EM(:,:,kk) = Psi_EM(:,:,kk+1) +v(:,:,kk)*(zf(1,1,kk+1)-zf(1,1,kk));
end

nW=1;
nE=nx;

for k=1:nz
for j=1:ny
thetay(j,k+1)=nanmean(squeeze(theta(nW:nE,j,k)));
%psiGM(j,k)=sum(Psiy_GM(nW:nE,j,k).*DXG(nW:nE,j));
psiEM(j,k)=sum(Psi_EM(nW:nE,j,k).*dxg(nW:nE,j));
end
end


