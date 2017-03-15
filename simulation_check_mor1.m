clear

flist = {'r056q','mor1'};
for ii=1:length(flist)
url='http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/aretxabaleta/floc/ocean_his_';
%url='http://geoport.whoi.edu/thredds/dodsC/usgs/data1/aretxabaleta/FLOC/gls_qcases/ocean_his_';
nc = ncgeodataset([url,flist{ii},'.nc'])
%%
run_name=flist{ii};
% read vertical grid parameters
Vtransform = nc{'Vtransform'}(:);
Vstretching = nc{'Vstretching'}(:);

s_rho = nc{'s_rho'}(:);
s_w = nc{'s_w'}(:);
Cs_r = nc{'Cs_r'}(:);
Cs_w = nc{'Cs_w'}(:);
N = length(s_rho);
Np = length(s_w);

theta_s = nc{'theta_s'}(:)
theta_b = nc{'theta_b'}(:)
depth_c = nc{'hc'}(:)

a=theta_s;
b=theta_b;
sr = s_rho;
C = (1-b)*sinh(a*sr)/sinh(a) + b*[tanh(a*(sr+0.5))/(2*tanh(0.5*a)) - 0.5];

%% read water depth
h = nc{'h'}(3,4);
hc = nc{'hc'}(:);
zeta = nc{'zeta'}(:,3,4);

z=squeeze(set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
                1, h, 0,0))';
zw=squeeze(set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
                 5, h, 0,0))';
             
time = nc{'ocean_time'}(:);
nt = length(time);
nz = length(z);
nzw = length(zw);
dz = diff(zw)

dz2d = repmat(dz,[nt,1]);
%%
fdiam = 1e6*nc{'Sd50'}(:)
ws = 1e3*nc{'Wsed'}(:)

% check mass conservation of suspended NCS classes only
% (last class in these runs is sand)
clear mud
ncs = length(ws)-1;
for n=1:ncs
   ncname = sprintf('mud_%02d',n)
   mud(n,:,:)=squeeze(nc{ncname}(:,:,3,4));
end
muds = squeeze(sum(mud));
mmud = muds.*dz2d;

initial_mass = sum(mmud(1,:));
final_mass = sum(mmud(nt-1,:));
summud_ts = sum(mmud,2);
max_mud_change = max( abs( summud_ts-summud_ts(1) ))

%%
% make 2D arrays of time and depths
t2d = repmat(time,[1,nz]);
tw2d = repmat(time,[1,nz+1]);
z2d = repmat(h+z,[nt,1]);
zw2d = repmat(h+zw,[nt,1]);

%%
% calculate size- and settling-velocity weighted floc dists
%eps = 1.e-8;
ws_av = squeeze(sum((mud+eps).*repmat(ws(1:ncs),[1,nt,nz]),1)./sum((mud+eps),1));
fdiam_av = squeeze(sum((mud+eps).*repmat(fdiam(1:ncs),[1,nt,nz]),1)./sum((mud+eps),1));

mconc = sum(mud,3);

%%
% plot time series of floc classes at various depths
for i=1:ncs
    if i<6;
        cl(i,:)=[0,0,(i)*50]/255;
    elseif i>=6&&i<10
        cl(i,:)=[(i-5)*50,0,(i-5)*50]/255;
    else
        cl(i,:)=[(i-10)*50,(i-10)*50,0]/255;
    end
    if ncs>16
        if i<7;
            cl(i,:)=[0,0,(i)*40]/255;
        elseif i>=7&&i<15
            cl(i,:)=[(i-6)*30,0,(i-6)*30]/255;
        else
            cl(i,:)=[(i-15)*30,(i-15)*30,0]/255;
        end
    end 
end

clear ho

%%Calculate and plot turbulence and mixing info
tke = nc{'tke'}(:,:,3,4);
gls = nc{'gls'}(:,:,3,4);
akv_bak = nc{'Akv_bak'}(:);
akt_bak = nc{'Akt_bak'}(:);
akv = nc{'AKv'}(:,:,3,4);
nueau =      1.5e-6;
gls_p =     -1.000;  %gls_p           GLS stability exponent.
gls_m =      0.500;  %gls_m           GLS turbulent kinetic energy exponent.
gls_n =     -1.000;  %gls_n           GLS turbulent length scale exponent.
gls_cmu0 =   5.4770e-01; %            GLS stability coefficient.
 
exp1 = 3.0+gls_p/gls_n;
exp2 = 1.5+gls_m/gls_n;
exp3 = -1.0/gls_n;
diss = gls_cmu0^exp1.*tke.^exp2.*gls.^exp3;
Gval=sqrt(diss/nueau)

u = nc{'u'}(:, :, 3, 4);

% read a 2D water level (height of ocean surface) at all time steps
%eta = nc.variables['zeta'][:, 3, 4]
% read bottom roughness zoN
zob = nc{'Zob'}(:);
zoa = nc{'Zo_app'}(:, 3, 4);
zoa(zoa<1e-6) = zob;
bustr = nc{'bustr'}(:, 3, 4);
bustrcwmax = nc{'bustrcwmax'}(:, 3, 4);
ubar = nc{'ubar'}(:, 3, 4);
%vbar = nc{'vbar'}(:, 3, 4);

%%

% Calculate G on rho 
Gc = 0.5*(Gval(:,1:end-1)+Gval(:,2:end));
if size(ws_av,1)>=10
    nts = 10;
else
    nts = size(ws_av,1)-1;
end


ws_av_mn = mean(ws_av(end-nts:end,:),1);
ws_av_mn2 = mean(ws_av,2);
ws_av_sd = std(ws_av(end-nts:end,:),1);
fdiam_av_mn = mean(fdiam_av(end-nts:end,:),1);
fdiam_av_sd = std(fdiam_av(end-nts:end,:),1);
tconc = squeeze(sum(mud,1));
tconc_mn = mean(tconc(end-nts:end,:),1);
tconc_sd = std(tconc(end-nts:end,:),1);
Gc_mn = mean(Gc(end-nts:end,:),1);
Gc_sd = std(Gc(end-nts:end,:),1);

nf = 2.;
m = .2;
q = (nf-1.)/(2.*m);

%%
figure(12);clf
set(gcf,'PaperPosition',[.5,.5,10,10]);wysiwyg
subplot (3,3,1)
pcolor(t2d/3600.,z2d,squeeze(sum(mud,1)));
shading flat
caxis([0,3])
title([' Concentration (kg/m^2) '])
colorbar;
set(gca,'FontSize',14,'FontWeight','bold')
subplot(3,3,4)
pcolor(t2d/3600.,z2d,fdiam_av/1000);
shading flat
caxis([4,4000]/1000)
ylabel('Elevation (m)')
title('Average diameter (mm) ')
colorbar;
set(gca,'FontSize',14,'FontWeight','bold')
subplot(3,3,7)
pcolor(t2d/3600.,z2d,ws_av);
shading flat;
caxis([0.1,10])
colorbar;
xlabel('Time (hrs)');
title('Settling velocity (mm/s) ')
set(gca,'FontSize',14,'FontWeight','bold')
subplot (3,3,2)
line(tconc_mn,h+z,'Color','k','LineWidth',1.5);
axis tight
title([' Concentration (kg/m^2) '])
set(gca,'FontSize',14,'FontWeight','bold')
subplot (3,3,5)
line(fdiam_av_mn/1000,h+z,'Color','k','LineWidth',1.5);
axis tight
title([' diameter (mm) '])
set(gca,'FontSize',14,'FontWeight','bold')
subplot (3,3,8)
line(ws_av_mn,h+z,'Color','k','LineWidth',1.5);
axis tight
title([' settling velocity (mm/s)  '])
set(gca,'FontSize',14,'FontWeight','bold')

subplot(1,3,3)
ap=get(gca,'Position');
for i=1:ncs
    ho(i)=semilogx(squeeze(mud(i,nt,:))+eps,h+z,'Color',cl(i,:),'LineWidth',1.5); hold on
end
ho(i+1)=semilogx(squeeze(sum(mud(:,nt-1,:),1)),h+z,'Color','k','LineStyle','--','LineWidth',2);
hl=legend(ho,[num2str(fdiam(1:ncs),'%1.0f');'Tota'],'Location','BestOutside');
set(hl,'FontSize',10,'FontWeight','bold','Position',[0.885,0.4,0.092,0.47]);
axis([0.0  20   0.01  12]);
set(gca,'Position',[0.67   0.11   0.213405797101449   0.815])
set(gca,'XTick',[1e-15,1e-12,1e-9,1e-6,1e-3,1e0],'FontSize',14,'FontWeight','bold')
xlabel('Concentration (kg/m3)')

eval(['print -dpng case_',run_name,'_summary.png'])


end
