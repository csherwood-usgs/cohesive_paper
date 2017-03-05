clear

flist = {'r056q','r057q','r058q','r059q','r060q','r061q',...
         'r062q','r063q','r064q','r065q','r066q','r067q','r068q','r069q','r070q','r071q',...
         'r072q','r073q','r074q','r075q','r076q','r077q'};
for ii=1:length(flist)
url='http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/aretxabaleta/floc/ocean_his_';
%url='http://geoport.whoi.edu/thredds/dodsC/usgs/data1/aretxabaleta/FLOC/gls_qcases/ocean_his_';
nc = ncgeodataset([url,flist{ii},'.nc'])
%%
run_name=flist{ii}(1:5);
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

if( max_mud_change > 1e-8 )
    % Make a plot if mass is not conserved5
    figure(1);clf 
    line(time/3600, sum(mmud,2),'Color', 'b');
    line(time/3600, sum(mmud,2),'Color', 'k','LineStyle','--');
    legend('Suspended','Total');
    ylabel('Mass (kg/m2)')
    xlabel('Time (hrs)')
end 
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
subplot (3,1,1)
pcolor(t2d/3600.,z2d,squeeze(sum(mud,1)));
shading flat
caxis([0.,0.02])
title([run_name,' Concentration (kg/m2)'])
colorbar;
subplot(3,1,2)
pcolor(t2d/3600.,z2d,fdiam_av);
shading flat
caxis([20,850])
ylabel('Elevation (m)')
title('Average diameter (mm)')
colorbar;
subplot(3,1,3)
pcolor(t2d/3600.,z2d,ws_av);
shading flat;
caxis([0.1,1])
colorbar;
xlabel('Time (hrs)');
title('Average settling velocity (mm/s)')

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
if ii~=10;
    legend(ho,[num2str(fdiam(1:ncs),'%1.0f');'Tota'],0);
else
    legend(ho,[num2str(fdiam(1:ncs),'%1.0f');'Tot'],0);
end
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

figure(8);clf
subplot(4,1,1)
pcolor(tw2d/3600.,zw2d,(akv+akv_bak));shading flat
title([run_name,' Eddy viscosity (m2/s)'])
colorbar;
subplot(4,1,2)
pcolor(tw2d/3600.,zw2d,Gval);shading flat
caxis([0,5])
colorbar;
title('Turbulence Shear Rate G (m/s)')
subplot(4,1,3)
pcolor(tw2d/3600.,zw2d,(tke));shading flat
title('Turbulence kinetic energy (m2/s2)')
colorbar
subplot(4,1,4)
contourf(t2d/3600.,z2d,u);
caxis([-.8,0])
colorbar;
xlabel('Time (hrs)')
title('Velocity (m/s)')
eval(['print -dpng case_',run_name,'mixing.png'])

%%
figure(9);clf
subplot(3,1,1)
plot(time/3600.,ubar);
title('Depth-Averaged Velocity')
ylabel('(m/s)')
subplot(3,1,2)
plot(time/3600.,bustr);hold on
plot(time/3600.,bustrcwmax,'Color','r')
title('Shear Stress')
ylabel('(Pa)')
subplot(3,1,3)
plot(time/3600.,log10(zoa))
title('Apparent Roughness')
ylabel('(m)')
xlabel('Time (hrs)')
eval(['print -dpng case_',run_name,'stress.png'])

% Calculate G on rho 
Gc = 0.5*(Gval(:,1:end-1)+Gval(:,2:end));
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

% report final values
display('Means for last 10 timesteps')
display('i   z    diam    ws     conc   Gc')
display([nz,   z2d(1,nz),  fdiam_av_mn(nz), ws_av_mn(nz), tconc_mn(nz), Gc_mn(nz)])
display([nz/2,   z2d(1,nz/2),  fdiam_av_mn(nz/2), ws_av_mn(nz/2), tconc_mn(nz/2), Gc_mn(nz/2)])
display([2,   z2d(1,2),  fdiam_av_mn(2), ws_av_mn(2), tconc_mn(2), Gc_mn(2)])
display([1,   z2d(1,1),  fdiam_av_mn(1), ws_av_mn(1), tconc_mn(1), Gc_mn(1)])
display('Standard deviations for last 10 timesteps')
display('i   z    diam    ws     conc   Gc')
display([nz,   z2d(1,nz),  fdiam_av_sd(nz), ws_av_sd(nz), tconc_sd(nz), Gc_sd(nz)])
display([nz/2,   z2d(1,nz/2),  fdiam_av_sd(nz/2), ws_av_sd(nz/2), tconc_sd(nz/2), Gc_sd(nz/2)])
display([2,   z2d(1,2),  fdiam_av_sd(2), ws_av_sd(2), tconc_sd(2), Gc_sd(2)])
display([1,   z2d(1,1),  fdiam_av_sd(1), ws_av_sd(1), tconc_sd(1), Gc_sd(1)])

nf = 2.;
m = .2;
q = (nf-1.)/(2.*m);
figure(9);clf
plot(tconc_mn./Gc_mn.^q,(1e-3*fdiam_av_mn).^(2*q),'.')
title(['Winterwerp plot: nf=',num2str(nf),' m=',num2str(m),' q=',num2str(q)])
ylabel('Diameter ^(2q)')
xlabel('Conc/G^q')
eval(['print -dpng case_',run_name,'winterwerp.png'])

figure(10);clf
title('Equilibrium Floc Size')
ylabel('Average Diameter (um)')
xlabel('Turbulence Shear Rate G (m/s)')
plot(Gc_mn,fdiam_av_mn,'.')
set(gca,'XLim',[0, 6])
eval(['print -dpng case_',run_name,'diam_v_G.png'])
eval(['save(''',run_name,'gc.txt'',' '''Gc_mn'',''-ascii'')'])
eval(['save(''',run_name,'diam.txt'',' '''fdiam_av_mn'',''-ascii'')'])
eval(['save(''',run_name,'ws.txt'',' '''ws_av_mn'',''-ascii'')'])
eval(['save(''',run_name,'mconc.txt'',' '''tconc_mn'',''-ascii'')'])
eval(['save (''',run_name,'_all.mat'',' '''Gc'',' '''fdiam_av'',' '''ws_av'',' '''tconc'',''-mat'')'])

end
%% plot results from several runs
% read data from GLS run, which has diffusivity for conc enabled
% G=zeros(length(flist),length(Gc_mn));
% d=zeros(size(G));
% c=zeros(size(G));
% ws=zeros(size(G));
for ii=1:length(flist)
    run_name=flist{ii}(1:5);
    G{ii} =load([run_name,'gc.txt']);
    d{ii}=load([run_name,'diam.txt']);
    wsa{ii}=load([run_name,'ws.txt']);
    c{ii}=load([run_name,'mconc.txt']);
end
for ii=1:length(flist)
    run_name=flist{ii}(1:5);
    load([run_name,'_all.mat'],'Gc');
    Gall{ii} =Gc(10:end,:);
    load([run_name,'_all.mat'],'fdiam_av');
    dall{ii}=fdiam_av(10:end,:);
    load([run_name,'_all.mat'],'ws_av');
    wsaall{ii}=ws_av(10:end,:);
    load([run_name,'_all.mat'],'tconc');
    call{ii}=tconc(10:end,:);
end

%%
% 56q 
%  settling and gls mixing
%  concentration 1 kg/m3 in class 3 (~0.04mm)
%  atm pressure grad 0.04 Pa
% 
% 57q 
%  settling and gls mixing
%  concentration 0.2 kg/m3 in class 3 (~0.04mm)
%  atm pressure grad 0.04 Pa
% 
% 58q 
%  settling and gls mixing
%  concentration 5 kg/m3 in class 3 (~0.04mm)
%  atm pressure grad 0.04 Pa
% 
% 59q 
%  settling and gls mixing
%  concentration 1 kg/m3 in class 3 (~0.04mm)
%  atm pressure grad 0.04 Pa
%  fractal dimension 1.9
%
% 60q 
%  settling and gls mixing
%  concentration 1 kg/m3 in class 3 (~0.04mm)
%  atm pressure grad 0.04 Pa
%  fractal dimension 2.1
% 
% 61q 
%  settling and gls mixing
%  concentration 1 kg/m3 in class 3 (~0.04mm)
%  atm pressure grad 0.04 Pa
%  fractal dimension 2.0
%  20 vertical layers
% 
% 62q 
%  settling and gls mixing
%  concentration 1 kg/m3 in class 3 (~0.04mm)
%  atm pressure grad 0.04 Pa
%  fractal dimension 2.0
%  100 vertical layers
% 
% 63q 
%  settling and gls mixing
%  concentration 1 kg/m3 in class 1 (0.02mm)
%  atm pressure grad 0.04 Pa
%  fractal dimension 2.0
%  50 vertical layers
% 
% 64q 
%  settling and gls mixing
%  concentration 1 kg/m3 in class 5 (~0.07mm)
%  atm pressure grad 0.04 Pa
%  fractal dimension 2.0
%  50 vertical layers
% 
% 65q 
%  settling and gls mixing
%  concentration 1 kg/m3 in class 2 (~0.04mm)
%  atm pressure grad 0.04 Pa
%  fractal dimension 2.0
%  50 vertical layers
%  7 cohesive classes
% 
% 66q 
%  settling and gls mixing
%  concentration 1 kg/m3 in class 6 (~0.04mm)
%  atm pressure grad 0.04 Pa
%  fractal dimension 2.0
%  50 vertical layers
%  22 cohesive classes


%% plot all results
figure(12);clf
%plt.title(r'Equilibrium Floc Size',fontsize=16)
ylabel('Average Diameter D (\mum)','FontSize',14)
xlabel('Turbulence Shear Rate G (s^{-1})','FontSize',14)
%plot(G[0,:],d[0,:],'.',markersize=10,label="1 kg/m^3")
line(G{3},d{3},'Marker','.','Color',[.5,0,0],'LineStyle','none')
line(G{1},d{1},'Marker','.','Color',[.8,.5,0],'LineStyle','none')
line(G{2},d{2},'Marker','.','Color',[.3,.3,.3],'LineStyle','none')
line(G{4},d{4},'Marker','.','Color',[.6,.6,.6],'LineStyle','none')
line(G{5}, d{5},'Marker','.','Color','g','LineStyle','none')
line(G{6}, d{6},'Marker','.','Color','b','LineStyle','none')
line(G{7}, d{7},'Marker','.','Color','k','LineStyle','none')
axis([0, 20,0,1500])
legend('5 kg/m^3','1 kg/m^3','0.2 kg/m^3','1 kg/m^3 nf=1.9','1 kg/m^3 nf=2.1',...
       '1 kg/m^3 N=20','1 kg/m^3 N=100' ,0)
print -dpng Multiple_cases_steady.png

%%

nf = 2.;
m = .2;
q = (nf-1.)/(2.*m);
q19 = (1.9-1.)/(2.*m);
q21 = (2.1-1.)/(2.*m);

figure(13);clf
%plot(tconc_mn/Gc_mn**q,(1e-3*fdiam_av_mn)**(2*q),'.')
%plot( c[0]/G[0]**q, (1e-3*d[0])**(2.*q),'.',label="1 kg/m^3")
line(c{3}./G{3}.^q, d{3}.^(2*q),'Marker','.','Color',[.5,0,0],'LineStyle','none')
line(c{1}./G{1}.^q, d{1}.^(2*q),'Marker','.','Color',[.8,.5,0],'LineStyle','none')
line(c{2}./G{2}.^q, d{2}.^(2*q),'Marker','.','Color',[.3,.3,.3],'LineStyle','none')
line(c{4}./G{4}.^q19, d{4}.^(2*q19),'Marker','.','Color',[.6,.6,.6],'LineStyle','none')
line(c{5}./G{5}.^q21, d{5}.^(2*q21),'Marker','.','Color','g','LineStyle','none')
line(c{6}./G{6}.^q, d{6}.^(2*q),'Marker','.','Color','b','LineStyle','none')
line(c{7}./G{7}.^q, d{7}.^(2*q),'Marker','.','Color','k','LineStyle','none')
title(['Winterwerp plot: nf=',num2str(nf),' m=',num2str(m),' q=',num2str(q)])
ylabel('Diameter ^(2q)')
xlabel('Conc/G^q')
legend('5 kg/m^3','1 kg/m^3','0.2 kg/m^3','1 kg/m^3 nf=1.9','1 kg/m^3 nf=2.1',...
       '1 kg/m^3 N=20','1 kg/m^3 N=100' ,0)
print -dpng Winterwerp_Multiple_cases_steady.png
%%
figure(23);clf
%plot(tconc_mn/Gc_mn**q,(1e-3*fdiam_av_mn)**(2*q),'.')
%plot( c[0]/G[0]**q, (1e-3*d[0])**(2.*q),'.',label="1 kg/m^3")
line(c{3}./G{3}.^q, d{3}.^(2*q),'Marker','.','Color',[.5,0,0],'LineStyle','none')
line(c{1}./G{1}.^q, d{1}.^(2*q),'Marker','.','Color',[.8,.5,0],'LineStyle','none')
line(c{2}./G{2}.^q, d{2}.^(2*q),'Marker','.','Color',[.3,.3,.3],'LineStyle','none')
line(c{4}./G{4}.^q19, d{4}.^(2*q19),'Marker','.','Color',[.6,.6,.6],'LineStyle','none')
line(c{5}./G{5}.^q21, d{5}.^(2*q21),'Marker','.','Color','g','LineStyle','none')
line(c{6}./G{6}.^q, d{6}.^(2*q),'Marker','.','Color','b','LineStyle','none')
line(c{7}./G{7}.^q, d{7}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{8}./G{8}.^q, d{8}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{9}./G{9}.^q, d{9}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{10}./G{10}.^q, d{10}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{11}./G{11}.^q, d{11}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{12}./G{12}.^q, d{12}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{13}./G{13}.^q, d{13}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{14}./G{14}.^q, d{14}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{15}./G{15}.^q, d{15}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{16}./G{16}.^q, d{16}.^(2*q),'Marker','.','Color','b','LineStyle','none')
line(c{17}./G{17}.^q, d{17}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{18}./G{18}.^q, d{18}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{19}./G{19}.^q, d{19}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{20}./G{20}.^q, d{20}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{21}./G{21}.^q, d{21}.^(2*q),'Marker','.','Color','k','LineStyle','none')
line(c{22}./G{22}.^q, d{22}.^(2*q),'Marker','.','Color','k','LineStyle','none')
title(['Winterwerp plot: nf=',num2str(nf),' m=',num2str(m),' q=',num2str(q)])
ylabel('Diameter ^(2q)')
xlabel('Conc/G^q')
legend('5 kg/m^3','1 kg/m^3','0.2 kg/m^3','1 kg/m^3 nf=1.9','1 kg/m^3 nf=2.1',...
       '1 kg/m^3 N=20','1 kg/m^3 N=100' ,0)
print -dpng Winterwerp_Multiple_cases_allo_steady.png
%%
figure(133);clf
%plot(tconc_mn/Gc_mn**q,(1e-3*fdiam_av_mn)**(2*q),'.')
%plot( c[0]/G[0]**q, (1e-3*d[0])**(2.*q),'.',label="1 kg/m^3")
a1=line(call{3}./Gall{3}.^q, dall{3}.^(2*q),'Marker','.','Color',[.5,0,0],'LineStyle','none');%58q
a2=line(call{1}./Gall{1}.^q, dall{1}.^(2*q),'Marker','.','Color',[.8,.5,0],'LineStyle','none');%56q
a3=line(call{2}./Gall{2}.^q, dall{2}.^(2*q),'Marker','.','Color',[.3,.3,.3],'LineStyle','none');%57q
a4=line(call{4}./Gall{4}.^q19, dall{4}.^(2*q19),'Marker','.','Color',[.6,.6,.6],'LineStyle','none');%59q
a5=line(call{5}./Gall{5}.^q21, dall{5}.^(2*q21),'Marker','.','Color','g','LineStyle','none');%60q
a6=line(call{6}./Gall{6}.^q, dall{6}.^(2*q),'Marker','.','Color','b','LineStyle','none');%61q
a7=line(call{7}./Gall{7}.^q, dall{7}.^(2*q),'Marker','.','Color','k','LineStyle','none');%62q
title(['Winterwerp plot: nf=',num2str(nf),' m=',num2str(m),' q=',num2str(q)])
ylabel('Diameter ^(2q)')
xlabel('Conc/G^q')
legend([a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1)],...
    '5 kg/m^3','1 kg/m^3','0.2 kg/m^3','1 kg/m^3 nf=1.9','1 kg/m^3 nf=2.1',...
       '1 kg/m^3 N=20','1 kg/m^3 N=100' ,0)
print -dpng Winterwerp_Multiple_cases_all_steady.png
%%
figure(131);clf
%plot(tconc_mn/Gc_mn**q,(1e-3*fdiam_av_mn)**(2*q),'.')
%plot( c[0]/G[0]**q, (1e-3*d[0])**(2.*q),'.',label="1 kg/m^3")
a1=line(call{3}./Gall{3}.^q, dall{3}.^(2*q),'Marker','.','Color',[.5,0,0],'LineStyle','none');%58q
a2=line(call{1}./Gall{1}.^q, dall{1}.^(2*q),'Marker','.','Color',[.8,.5,0],'LineStyle','none');%56q
a3=line(call{2}./Gall{2}.^q, dall{2}.^(2*q),'Marker','.','Color',[.3,.3,.3],'LineStyle','none');%57q
%a4=line(call{4}./Gall{4}.^q19, dall{4}.^(2*q19),'Marker','.','Color',[.6,.6,.6],'LineStyle','none');%59q
%a5=line(call{5}./Gall{5}.^q21, dall{5}.^(2*q21),'Marker','.','Color','g','LineStyle','none');%60q
a6=line(call{6}./Gall{6}.^q, dall{6}.^(2*q),'Marker','.','Color','b','LineStyle','none');%61q
a7=line(call{7}./Gall{7}.^q, dall{7}.^(2*q),'Marker','.','Color','k','LineStyle','none');%62q
title(['Winterwerp plot: nf=',num2str(nf),' m=',num2str(m),' q=',num2str(q)])
ylabel('Diameter ^(2q)')
xlabel('Conc/G^q')
legend([a1(1),a2(1),a3(1),a6(1),a7(1)],...
    '5 kg/m^3','1 kg/m^3','0.2 kg/m^3',...
       '1 kg/m^3 N=20','1 kg/m^3 N=100' ,0)
print -dpng Winterwerp_Multiple_cases_all_nf2_steady.png
%%
figure(14);clf
xlabel('Average concentration ','FontSize',14)
ylabel('settling velocity (s^{-1})','FontSize',14)
line(c{3},wsa{3},'Marker','.','Color',[.8,.5,0],'LineStyle','none','MarkerSize',16)
line(c{1},wsa{1},'Marker','.','Color',[.3,.3,.3],'LineStyle','none','MarkerSize',16)
line(c{2},wsa{2},'Marker','.','Color',[.5,0,0],'LineStyle','none','MarkerSize',16)
line(c{4},wsa{4},'Marker','.','Color',[.6,.6,.6],'LineStyle','none','MarkerSize',16)
line(c{5},wsa{5},'Marker','.','Color','g','LineStyle','none','MarkerSize',16)
line(c{6},wsa{6},'Marker','.','Color','b','LineStyle','none','MarkerSize',16)
line(c{7},wsa{7},'Marker','.','Color','k','LineStyle','none','MarkerSize',16)
set(gca,'YScale','log','XScale','log');grid
%axis([0, 6,100,700])
legend('5 kg/m^3','1 kg/m^3','0.2 kg/m^3','1 kg/m^3 nf=1.9','1 kg/m^3 nf=2.1',...
       '1 kg/m^3 N=20','1 kg/m^3 N=100' ,0)
print -dpng Multiple_cases_diam_WS_steady.png
%%
figure(15);clf
set(gcf,'PaperPosition',[.5,.5,10,10]);wysiwyg
xlabel('Average concentration ','FontSize',14)
ylabel('settling velocity (s^{-1})','FontSize',14)
set(gca,'YScale','log','XScale','log');grid
a1=line(call{3}',wsaall{3}','Marker','.','Color',[.8,.5,0],'LineStyle','none','MarkerSize',16); %58q
a2=line(call{1}',wsaall{1}','Marker','.','Color',[.3,.3,.3],'LineStyle','none','MarkerSize',16); %56q
a3=line(call{2}',wsaall{2}','Marker','.','Color',[.5,0,0],'LineStyle','none','MarkerSize',16);   %57q
a4=line(call{4}',wsaall{4}','Marker','.','Color',[.6,.6,.6],'LineStyle','none','MarkerSize',16); %59q
a5=line(call{5}',wsaall{5}','Marker','.','Color','g','LineStyle','none','MarkerSize',16);        %60q
a6=line(call{6}',wsaall{6}','Marker','.','Color','b','LineStyle','none','MarkerSize',16);        %61q
a7=line(call{7}',wsaall{7}','Marker','.','Color','k','LineStyle','none','MarkerSize',16);        %62q
a8=line(call{8}',wsaall{8}','Marker','.','Color',[0,.5,0],'LineStyle','none','MarkerSize',16);   %63q
a9=line(call{9}',wsaall{9}','Marker','.','Color',[0,.5,.5],'LineStyle','none','MarkerSize',16);  %64q
a10=line(call{10}',wsaall{10}','Marker','.','Color','c','LineStyle','none','MarkerSize',16);     %65q
a11=line(call{11}',wsaall{11}','Marker','.','Color','y','LineStyle','none','MarkerSize',16);     %66q
%axis([1e-1, 10,0.07,0.3])
set(gca,'YTick',[0.08,0.09,0.1,0.2,0.3])
legend([a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1),a9(1),a10(1),a11(1)],...
    '5 kg/m^3','1 kg/m^3','0.2 kg/m^3','1 kg/m^3 nf=1.9','1 kg/m^3 nf=2.1',...
       '1 kg/m^3 N=20','1 kg/m^3 N=100','1 kg/m^3 cl_o=1','1 kg/m^3 cl_o=5','1 kg/m^3 NCS=7','1 kg/m^3 NCS=22',...
    0);
print -dpng Multiple_cases_diam_WS_all_steady.png
%%
figure(16);clf
set(gcf,'PaperPosition',[.5,.5,10,10]);wysiwyg
xlabel('Average concentration ','FontSize',14)
ylabel('settling velocity (s^{-1})','FontSize',14)
set(gca,'YScale','log','XScale','log');grid
%a1=line(call{3}',wsaall{3}','Marker','.','Color',[.8,.5,0],'LineStyle','none','MarkerSize',16); %58q
a2=line(call{1}',wsaall{1}','Marker','.','Color',[.3,.3,.3],'LineStyle','none','MarkerSize',16); %56q
%a3=line(call{2}',wsaall{2}','Marker','.','Color',[.5,0,0],'LineStyle','none','MarkerSize',16);   %57q
a4=line(call{4}',wsaall{4}','Marker','.','Color',[.6,.6,.6],'LineStyle','none','MarkerSize',16); %59q
a5=line(call{5}',wsaall{5}','Marker','.','Color','g','LineStyle','none','MarkerSize',16);        %60q
a6=line(call{6}',wsaall{6}','Marker','.','Color','b','LineStyle','none','MarkerSize',16);        %61q
a7=line(call{7}',wsaall{7}','Marker','.','Color','k','LineStyle','none','MarkerSize',16);        %62q
a8=line(call{8}',wsaall{8}','Marker','.','Color',[0,.5,0],'LineStyle','none','MarkerSize',16);   %63q
a9=line(call{9}',wsaall{9}','Marker','.','Color',[0,.5,.5],'LineStyle','none','MarkerSize',16);  %64q
a10=line(call{10}',wsaall{10}','Marker','.','Color','c','LineStyle','none','MarkerSize',16);     %65q
a11=line(call{11}',wsaall{11}','Marker','.','Color','y','LineStyle','none','MarkerSize',16);     %66q
a10=line(call{10}',wsaall{10}','Marker','.','Color','c','LineStyle','none','MarkerSize',16);     %65q
%axis([6e-1, 2,0.07,0.3])
set(gca,'YTick',[0.08,0.09,0.1,0.2,0.3])
legend([a2(1),a4(1),a5(1),a6(1),a7(1),a8(1),a9(1),a10(1),a11(1)],...
    '1 kg/m^3','1 kg/m^3 nf=1.9','1 kg/m^3 nf=2.1',...
       '1 kg/m^3 N=20','1 kg/m^3 N=100','1 kg/m^3 cl_o=1','1 kg/m^3 cl_o=5','1 kg/m^3 NCS=7','1 kg/m^3 NCS=22',...
    0);
print -dpng Multiple_cases_diam_WS_1kg_steady.png
%%
figure(17);clf
set(gcf,'PaperPosition',[.5,.5,10,10]);wysiwyg
xlabel('Average concentration ','FontSize',14)
ylabel('settling velocity (mm s^{-1})','FontSize',14)
set(gca,'YScale','log','XScale','log');grid
a3=line(call{3}',wsaall{3}','Marker','.','Color',[.8,.5,0],'LineStyle','none','MarkerSize',16); %58q
a1=line(call{1}',wsaall{1}','Marker','.','Color',[.3,.3,.3],'LineStyle','none','MarkerSize',16);%56q
a2=line(call{2}',wsaall{2}','Marker','.','Color',[.5,0,0],'LineStyle','none','MarkerSize',16); %57q
a12=line(call{12}',wsaall{12}','Marker','.','Color',[0,0.5,0],'LineStyle','none','MarkerSize',16);
a13=line(call{13}',wsaall{13}','Marker','.','Color',[0,0.5,.5],'LineStyle','none','MarkerSize',16);
a14=line(call{14}',wsaall{14}','Marker','.','Color',[.5,0.5,0],'LineStyle','none','MarkerSize',16);
a15=line(call{15}',wsaall{15}','Marker','.','Color',[.8,.8,.8],'LineStyle','none','MarkerSize',16);
a16=line(call{16}',wsaall{16}','Marker','.','Color',[.5,.5,.5],'LineStyle','none','MarkerSize',16);
a17=line(call{17}',wsaall{17}','Marker','.','Color','r','LineStyle','none','MarkerSize',16);
a18=line(call{18}',wsaall{18}','Marker','.','Color','b','LineStyle','none','MarkerSize',16);
a19=line(call{19}',wsaall{19}','Marker','.','Color','k','LineStyle','none','MarkerSize',16);
a20=line(call{20}',wsaall{20}','Marker','.','Color','g','LineStyle','none','MarkerSize',16);
a21=line(call{21}',wsaall{21}','Marker','.','Color','c','LineStyle','none','MarkerSize',16);
a22=line(call{22}',wsaall{22}','Marker','.','Color','y','LineStyle','none','MarkerSize',16);
%axis([7e-2, 130,0.08,0.15])
set(gca,'YTick',[0.08,0.09,0.1,0.12,0.15])
legend([a12(1),a2(1),a13(1),a14(1),a15(1),a1(1),a16(1),a17(1),a3(1),a18(1),a19(1),a20(1),a21(1),a22(1)],...
    '0.1 kg/m^3','0.2 kg/m^3','0.3 kg/m^3','0.5 kg/m^3','0.7 kg/m^3','1 kg/m^3','2 kg/m^3',...
    '3 kg/m^3','5 kg/m^3','7 kg/m^3','10 kg/m^3','20 kg/m^3','50 kg/m^3','100 kg/m^3',...
    'Location','BestOutside');
print -dpng Multiple_cases_diam_WS_diffConc_steady.png


%%
figure(144);clf
xlabel('Average diameter ','FontSize',14)
ylabel('settling velocity (mm s^{-1})','FontSize',14)
a1=line(dall{3}, wsaall{3},'Marker','.','Color',[.5,0,0],'LineStyle','none');
a2=line(dall{1}, wsaall{1},'Marker','.','Color',[.8,.5,0],'LineStyle','none');
a3=line(dall{2}, wsaall{2},'Marker','.','Color',[.3,.3,.3],'LineStyle','none');
a4=line(dall{4}, wsaall{4},'Marker','.','Color',[.6,.6,.6],'LineStyle','none');
a5=line(dall{5}, wsaall{5},'Marker','.','Color','g','LineStyle','none');
a6=line(dall{6}, wsaall{6},'Marker','.','Color','b','LineStyle','none');
a7=line(dall{7}, wsaall{7},'Marker','.','Color','k','LineStyle','none');
a8=line(dall{8}, wsaall{8},'Marker','.','Color','c','LineStyle','none');
a9=line(dall{9}, wsaall{9},'Marker','.','Color','y','LineStyle','none');
a10=line(dall{10},wsaall{10},'Marker','.','Color','g','LineStyle','none');
a11=line(dall{11},wsaall{11},'Marker','.','Color','k','LineStyle','none');
a12=line(dall{12},wsaall{12},'Marker','.','Color','k','LineStyle','none');
a13=line(dall{13},wsaall{13},'Marker','.','Color','k','LineStyle','none');
a14=line(dall{14},wsaall{14},'Marker','.','Color','k','LineStyle','none');
a15=line(dall{15},wsaall{15},'Marker','.','Color','k','LineStyle','none');
a16=line(dall{16},wsaall{16},'Marker','.','Color','b','LineStyle','none');
a17=line(dall{17},wsaall{17},'Marker','.','Color','k','LineStyle','none');
a18=line(dall{18},wsaall{18},'Marker','.','Color','k','LineStyle','none');
a19=line(dall{19},wsaall{19},'Marker','.','Color','k','LineStyle','none');
a20=line(dall{20},wsaall{20},'Marker','.','Color',[127/255 212/255, 1],'LineStyle','none');
a21=line(dall{21},wsaall{21},'Marker','.','Color',[1,127/255 212/255],'LineStyle','none');
a22=line(dall{22},wsaall{22},'Marker','.','Color',[127/255 1 212/255],'LineStyle','none');%axtt
set(gca,'YScale','log','XScale','log');grid
%axis([0, 6,100,700])
legend([a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1),a9(1),a10(1),a11(1),a12(1),...
        a13(1),a14(1),a15(1),a16(1),a17(1),a18(1),a19(1),a20(1),a21(1),a22(1),],'5 kg/m^3','1 kg/m^3','0.2 kg/m^3','1 kg/m^3 nf=1.9','1 kg/m^3 nf=2.1',...
       '1 kg/m^3 N=20','1 kg/m^3 N=100' ,0)
print -dpng Multiple_cases_diam_WS_steady.png

