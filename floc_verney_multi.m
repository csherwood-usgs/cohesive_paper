clear

flist = {'verneyc1','verneyc3','verneyc5','verneyc7','verneyc9','verneyc11',...
    'verneyc13','verneyc15','verneyc17','verneyc19','verneyc21'};
for ii=1:length(flist)
url='http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/aretxabaleta/floc/ocean_his_';
nc = ncgeodataset([url,flist{ii},'.nc'])

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

figure(3);clf
pcolor(t2d/3600.,z2d,fdiam_av);
shading flat
caxis([20,1500])
ylabel('Elevation (m)')
title('Average diameter (mm)')
colorbar;
xlabel('Time (hrs)');
%title('Average settling velocity (mm/s)')

%%
display('i   z    diam    ws     conc')
display([nz,   z2d(1,nz),  fdiam_av(end,nz), ws_av(end,nz),sum(mud(:,nt,nz),1 )])
display([nz/2,   z2d(1,nz/2),  fdiam_av(end,nz/2), ws_av(end,nz/2),sum(mud(:,nt,nz/2),1 )])
display([2,   z2d(1,2),  fdiam_av(end,2), ws_av(end,2),sum(mud(:,nt,2),1 )])
display([1,   z2d(1,1),  fdiam_av(end,1), ws_av(end,1),sum(mud(:,nt,1),1 )])

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
        if i<6;
            cl(i,:)=[0,0,(i)*20]/255;
        elseif i>=6&&i<10
            cl(i,:)=[(i-5)*20,0,(i-5)*20]/255;
        else
            cl(i,:)=[(i-10)*20,(i-10)*20,0]/255;
        end
    end 
end
figure(6);clf
set(gcf,'PaperPosition',[.5,.5,8,10.5]);wysiwyg
subplot(4,1,1)
j=nz-1
lb = h+z(j);
for i=1:ncs
    m = mud(i,:,j);
    line(time/3600,m,'Color',cl(i,:));
	title([run_name,' ',num2str(lb,2),' mab'])
end    
subplot(4,1,2)
j=nz/2
lb = h+z(j);
for i=1:ncs
    m = mud(i,:,j);
    line(time/3600,m,'Color',cl(i,:));
	title([run_name,' ',num2str(lb,2),' mab'])
end    
subplot(4,1,3)
j=2
lb = h+z(j);
for i=1:ncs
    m = mud(i,:,j);
    line(time/3600,m,'Color',cl(i,:));
	title([run_name,' ',num2str(lb,2),' mab'])
end    
subplot(4,1,4)
for i=1:ncs
    m = sum(mud(i,:,:),3);
    ho(i)=line(time/3600,m,'Color',cl(i,:));
    title('Total')
end
legend(ho,num2str(fdiam,'%1.0f'),0)
eval(['print -dpng case_',run_name,'.png'])

clear ho
%%Final profiles
figure(7);clf
% plot individual size classes
for i=1:ncs
    ho(i)=semilogx(squeeze(mud(i,nt,:))+eps,h+z,'Color',cl(i,:)); hold on
end
ho(i+1)=semilogx(squeeze(sum(mud(:,nt-1,:),1)),h+z,'Color','k','LineStyle','--','LineWidth',2);
%if ii~=10;
    legend(ho,[num2str(fdiam(1:ncs),'%1.0f');'Tota'],0);
%else
%    legend(ho,[num2str(fdiam(1:ncs),'%1.0f');'Tot'],0);
%end
title([run_name,' Final Profiles'])
ylabel('Elevation (m)')
xlabel('Concentration (kg/m3)')
eval(['print -dpng case_',run_name,'final_profiles.png'])

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

%tidx = -1       % just get the final time step, for now.
% read a 3D temperature field at specified time step
%temp = nc.variables['temp'][tidx, :, :, :]
u = nc{'u'}(:, :, 3, 4);
%v = nc.variables['v'][tidx, :, :, :]

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

Gc=0.0*Gc;
jj=find(t2d<7201);
Gc(jj)=1;
jj=find(t2d>=7201&t2d<8401);
Gc(jj)=2;
jj=find(t2d>=8401&t2d<9601);
Gc(jj)=3;
jj=find(t2d>=9601&t2d<10801);
Gc(jj)=4;
jj=find(t2d>=10801&t2d<12601);
Gc(jj)=12;
jj=find(t2d>=12601&t2d<13801);
Gc(jj)=4;
jj=find(t2d>=13801&t2d<15001);
Gc(jj)=3;
jj=find(t2d>=15001&t2d<16201);
Gc(jj)=2;
jj=find(t2d>=16201&t2d<21601);
Gc(jj)=1;
jj=find(t2d>=21601&t2d<25201);
Gc(jj)=0;
jj=find(t2d>=25201&t2d<30601);
Gc(jj)=1;
jj=find(t2d>=30601&t2d<31801);
Gc(jj)=2;
jj=find(t2d>=31801&t2d<33001);
Gc(jj)=3;
jj=find(t2d>=33001&t2d<34201);
Gc(jj)=4;
jj=find(t2d>=34201&t2d<36001);
Gc(jj)=12;
jj=find(t2d>=36001&t2d<37201);
Gc(jj)=4;
jj=find(t2d>=37201&t2d<38401);
Gc(jj)=3;
jj=find(t2d>=38401&t2d<39601);
Gc(jj)=2;
jj=find(t2d>=39601&t2d<45001);
Gc(jj)=1;
jj=find(t2d>=45001&t2d<48601);
Gc(jj)=0;
jj=find(t2d>=48601&t2d<54001);
Gc(jj)=1;
jj=find(t2d>=54001&t2d<55201);
Gc(jj)=2;
jj=find(t2d>=55201&t2d<56401);
Gc(jj)=3;
jj=find(t2d>=56401&t2d<57601);
Gc(jj)=4;
jj=find(t2d>=57601);
Gc(jj)=12;


nts = 10;

ws_av_mn = mean(ws_av(end-10:end,:),1);
ws_av_mn2 = mean(ws_av,2);
ws_av_sd = std(ws_av(end-10:end,:),1);
fdiam_av_mn = mean(fdiam_av(end-10:end,:),1);
fdiam_av_sd = std(fdiam_av(end-10:end,:),1);
tconc = squeeze(sum(mud,1));
tconc_mn = mean(tconc(end-10:end,:),1);
tconc_sd = std(tconc(end-10:end,:),1);
Gc_mn = mean(Gc(end-10:end,:),1);
Gc_sd = std(Gc(end-10:end,:),1);

%%
figure(4);clf
subplot (2,1,1)
pcolor(t2d/60.,z2d,squeeze(Gc));
shading flat
caxis([0.,12])
title([' Turbulence Shear Rate G (s ^{-1})'])
colorbar;
subplot(2,1,2)
pcolor(t2d/60,z2d,fdiam_av);
shading flat
caxis([20,1500])
ylabel('Elevation (m)')
title('Average diameter (mm) ')
colorbar;
xlabel('Time (min)');
%%
figure(5);clf
subplot (2,1,1)
line(t2d(:,25)/60,squeeze(Gc(:,25)));
title([' Turbulence Shear Rate G (s ^{-1})'])
subplot(2,1,2)
line(t2d(:,25)/60,fdiam_av(:,50),'Color','r');
line(t2d(:,25)/60,fdiam_av(:,1),'Color','g');
line(t2d(:,25)/60,fdiam_av(:,25));
ylabel('Average diameter (mm) ')
xlabel('Time (min)');
eval(['print -dpng case_',run_name,'_G_D.png'])
%%
figure(8);clf
subplot (2,1,1)
line(t2d(:,25)/60,squeeze(Gc(:,25)));
title([' Turbulence Shear Rate G (s ^{-1})'])
subplot(2,1,2)
line(t2d(:,25)/60,squeeze(mud(:,:,25)));
ylabel('Average Concentration (kg m^{-3}) ')
xlabel('Time (min)');

%%
allr{ii}.tt=t2d(:,25)/60;
allr{ii}.fdiam_25=fdiam_av(:,25);
allr{ii}.G=squeeze(Gc(:,25));
allr{ii}.mud=squeeze(mud(:,:,25));

end

%%
lineStyles=linspecer(length(allr));

figure(55);clf
set(gcf,'PaperPosition',[.5,.5,10.5,8]);wysiwyg
% ha = tight_subplot(2,1,[.05 .03],[.1 .1],[.1 .18]);
% axes(ha(1));
subplot('Position',[0.1   0.675   0.72   0.225])
line(t2d(:,25)/60,squeeze(Gc(:,25)),'LineWidth',1.5);grid
title([' Turbulence Shear Rate G (s ^{-1})'])
%axes(ha(2));
set(gca,'XLim',[0,1060])
subplot('Position',[0.1   0.1   0.72   0.525])
for ii=1:length(allr)
    h1(ii)=line(allr{ii}.tt,allr{ii}.fdiam_25,'Color',lineStyles(ii,:),'LineWidth',1.5);
end
aa=get(gca,'Position');
set(gca,'YLim',[0,700])
ylabel('Average diameter (mm) ')
xlabel('Time (min)');grid
set(gca,'XLim',[0,1060])
hle=legend([h1],'Class 1','Class 3','Class 5','Class 7','Class 9','Class 11',...
    'Class 13','Class 15','Class 17','Class 19','Class 21',...
       'Location','BestOutside');
set(gca,'Position',aa)
eval(['print -dpng all_verney_G_D.png'])
%%
figure(56);clf
set(gcf,'PaperPosition',[.5,.5,10.5,8]);wysiwyg
subplot('Position',[0.1   0.675   0.72   0.225])
line(t2d(:,25)/60,squeeze(Gc(:,25)),'LineWidth',1.5);grid
title([' Turbulence Shear Rate G (s ^{-1})'])
set(gca,'XLim',[0,1060])
subplot('Position',[0.1   0.1   0.72   0.525])
for ii=1:length(allr)
    h1(ii)=line(allr{ii}.tt,allr{ii}.fdiam_25,'Color',lineStyles(ii,:),'LineWidth',1.5);
end
aa=get(gca,'Position');
set(gca,'YLim',[0,5000],'XLim',[0,1060])
ylabel('Average diameter (mm) ')
xlabel('Time (min)');grid
hle=legend([h1],'Class 1','Class 3','Class 5','Class 7','Class 9','Class 11',...
    'Class 13','Class 15','Class 17','Class 19','Class 21',...
       'Location','BestOutside');
set(gca,'Position',aa,'YScale','log')
eval(['print -dpng all_verney_G_D_log.png'])

