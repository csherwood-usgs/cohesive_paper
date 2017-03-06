clear

flist = {'noset056q','noset057q','noset058q','noset067q','noset068q','noset069q','noset070q','noset071q',...
         'noset072q','noset073q','noset074q','noset075q','noset076q','noset077q'};
for ii=1:length(flist)
url='http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/aretxabaleta/floc/ocean_his_';
%url='http://geoport.whoi.edu/thredds/dodsC/usgs/data1/aretxabaleta/FLOC/gls_qcases/ocean_his_';
nc = ncgeodataset([url,flist{ii},'.nc'])
%%
run_name=flist{ii}(1:9);
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
fdiam = 1e6*nc{'Sd50'}(1)
ws = 1e3*nc{'Wsed'}(1)

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
legend(ho,num2str(fdiam,'%1.0f'))
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
    legend(ho,[num2str(fdiam(1:ncs),'%1.0f');'Tota']);
%else
%    legend(ho,[num2str(fdiam(1:ncs),'%1.0f');'Tot'],0);
%end
title([run_name,' Final Profiles'])
ylabel('Elevation (m)')
xlabel('Concentration (kg/m3)')
eval(['print -dpng case_',run_name,'final_profiles.png'])

%% Calculate and plot turbulence and mixing info
tke = nc{'tke'}(:,:,3,4);
gls = nc{'gls'}(:,:,3,4);
akv_bak = nc{'Akv_bak'}(1);
akt_bak = nc{'Akt_bak'}(1);
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

%% read a 2D water level (height of ocean surface) at all time steps
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
    run_name=flist{ii}(1:9);
    G{ii} =load([run_name,'gc.txt']);
    d{ii}=load([run_name,'diam.txt']);
    wsa{ii}=load([run_name,'ws.txt']);
    c{ii}=load([run_name,'mconc.txt']);
end
for ii=1:length(flist)
    run_name=flist{ii}(1:9);
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
figure(122);clf
ylabel('Average Diameter D (\mum)','FontSize',14)
xlabel('Turbulence Shear Rate G (s^{-1})','FontSize',14)
a14=line(G{14}, d{14},'Marker','.','MarkerSize',16,'Color',[.5,0,.5],'LineStyle','none');%100 %77
a13=line(G{13}, d{13},'Marker','.','MarkerSize',16,'Color',[0,0,.5],'LineStyle','none');%50   %76
a12=line(G{12}, d{12},'Marker','.','MarkerSize',16,'Color',[0,.5,.5],'LineStyle','none');%20  %75
a11=line(G{11}, d{11},'Marker','.','MarkerSize',16,'Color',[.8,.8,.8],'LineStyle','none');%10 %74
a10=line(G{10}, d{10},'Marker','.','MarkerSize',16,'Color','y','LineStyle','none');%7         %73
a3=line(G{3}, d{3},'Marker','.','MarkerSize',16,'Color',[.5,0,0],'LineStyle','none');%5       %58
a9=line(G{9}, d{9},'Marker','.','MarkerSize',16,'Color','r','LineStyle','none');%3            %72
a8=line(G{8}, d{8},'Marker','.','MarkerSize',16,'Color','c','LineStyle','none');%2            %71
a1=line(G{1}, d{1},'Marker','.','MarkerSize',16,'Color',[.8,.5,0],'LineStyle','none');%1      %56
a7=line(G{7}, d{7},'Marker','.','MarkerSize',16,'Color','k','LineStyle','none');%0.7          %70
a6=line(G{6}, d{6},'Marker','.','MarkerSize',16,'Color','b','LineStyle','none');%0.5          %69
a5=line(G{5}, d{5},'Marker','.','MarkerSize',16,'Color','g','LineStyle','none');%0.3          %68
a2=line(G{2}, d{2},'Marker','.','MarkerSize',16,'Color',[.3,.3,.3],'LineStyle','none');%0.2   %57
a4=line(G{4}, d{4},'Marker','.','MarkerSize',16,'Color',[.6,.6,.6],'LineStyle','none');%0.1   %67
axis([0, 10,0,1500])
legend([a4(1),a2(1),a5(1),a6(1),a7(1),a1(1),a8(1),a9(1),a3(1),a10(1),a11(1),a12(1),a13(1),a14(1)],...
       '0.1 kg/m^3','0.2 kg/m^3','0.3 kg/m^3','0.5 kg/m^3','0.7 kg/m^3','1 kg/m^3',...
       '2 kg/m^3','3 kg/m^3','5 kg/m^3','7 kg/m^3','10 kg/m^3','20 kg/m^3','50 kg/m^3','100 kg/m^3' ,...
       'Location','BestOutside')
print -dpng Multiple_cases_noset_G_D.png
%%
figure(123);clf
ylabel('Average Concentration (kg m^{-3})','FontSize',14)
xlabel('Turbulence Shear Rate G (s^{-1})','FontSize',14)
a14=line(G{14}, c{14},'Marker','.','MarkerSize',16,'Color',[.5,0,.5],'LineStyle','none');%100 %77
a13=line(G{13}, c{13},'Marker','.','MarkerSize',16,'Color',[0,0,.5],'LineStyle','none');%50   %76
a12=line(G{12}, c{12},'Marker','.','MarkerSize',16,'Color',[0,.5,.5],'LineStyle','none');%20  %75
a11=line(G{11}, c{11},'Marker','.','MarkerSize',16,'Color',[.8,.8,.8],'LineStyle','none');%10 %74
a10=line(G{10}, c{10},'Marker','.','MarkerSize',16,'Color','y','LineStyle','none');%7         %73
a3=line(G{3}, c{3},'Marker','.','MarkerSize',16,'Color',[.5,0,0],'LineStyle','none');%5       %58
a9=line(G{9}, c{9},'Marker','.','MarkerSize',16,'Color','r','LineStyle','none');%3            %72
a8=line(G{8}, c{8},'Marker','.','MarkerSize',16,'Color','c','LineStyle','none');%2            %71
a1=line(G{1}, c{1},'Marker','.','MarkerSize',16,'Color',[.8,.5,0],'LineStyle','none');%1      %56
a7=line(G{7}, c{7},'Marker','.','MarkerSize',16,'Color','k','LineStyle','none');%0.7          %70
a6=line(G{6}, c{6},'Marker','.','MarkerSize',16,'Color','b','LineStyle','none');%0.5          %69
a5=line(G{5}, c{5},'Marker','.','MarkerSize',16,'Color','g','LineStyle','none');%0.3          %68
a2=line(G{2}, c{2},'Marker','.','MarkerSize',16,'Color',[.3,.3,.3],'LineStyle','none');%0.2   %57
a4=line(G{4}, c{4},'Marker','.','MarkerSize',16,'Color',[.6,.6,.6],'LineStyle','none');%0.1   %67
axis([0, 10,0,100])
legend([a4(1),a2(1),a5(1),a6(1),a7(1),a1(1),a8(1),a9(1),a3(1),a10(1),a11(1),a12(1),a13(1),a14(1)],...
       '0.1 kg/m^3','0.2 kg/m^3','0.3 kg/m^3','0.5 kg/m^3','0.7 kg/m^3','1 kg/m^3',...
       '2 kg/m^3','3 kg/m^3','5 kg/m^3','7 kg/m^3','10 kg/m^3','20 kg/m^3','50 kg/m^3','100 kg/m^3' ,...
       'Location','BestOutside')
print -dpng Multiple_cases_noset_G_C.png
%%
figure(124);clf
ylabel('Average Diameter D (\mum)','FontSize',14)
xlabel('Average Concentration (kg m^{-3})','FontSize',14)
a14=line(c{14}, d{14},'Marker','.','MarkerSize',16,'Color',[.5,0,.5],'LineStyle','none');%100 %77
a13=line(c{13}, d{13},'Marker','.','MarkerSize',16,'Color',[0,0,.5],'LineStyle','none');%50   %76
a12=line(c{12}, d{12},'Marker','.','MarkerSize',16,'Color',[0,.5,.5],'LineStyle','none');%20  %75
a11=line(c{11}, d{11},'Marker','.','MarkerSize',16,'Color',[.8,.8,.8],'LineStyle','none');%10 %74
a10=line(c{10}, d{10},'Marker','.','MarkerSize',16,'Color','y','LineStyle','none');%7         %73
a3=line(c{3}, d{3},'Marker','.','MarkerSize',16,'Color',[.5,0,0],'LineStyle','none');%5       %58
a9=line(c{9}, d{9},'Marker','.','MarkerSize',16,'Color','r','LineStyle','none');%3            %72
a8=line(c{8}, d{8},'Marker','.','MarkerSize',16,'Color','c','LineStyle','none');%2            %71
a1=line(c{1}, d{1},'Marker','.','MarkerSize',16,'Color',[.8,.5,0],'LineStyle','none');%1      %56
a7=line(c{7}, d{7},'Marker','.','MarkerSize',16,'Color','k','LineStyle','none');%0.7          %70
a6=line(c{6}, d{6},'Marker','.','MarkerSize',16,'Color','b','LineStyle','none');%0.5          %69
a5=line(c{5}, d{5},'Marker','.','MarkerSize',16,'Color','g','LineStyle','none');%0.3          %68
a2=line(c{2}, d{2},'Marker','.','MarkerSize',16,'Color',[.3,.3,.3],'LineStyle','none');%0.2   %57
a4=line(c{4}, d{4},'Marker','.','MarkerSize',16,'Color',[.6,.6,.6],'LineStyle','none');%0.1   %67
axis([0, 10,0,1500])
legend([a4(1),a2(1),a5(1),a6(1),a7(1),a1(1),a8(1),a9(1),a3(1),a10(1),a11(1),a12(1),a13(1),a14(1)],...
       '0.1 kg/m^3','0.2 kg/m^3','0.3 kg/m^3','0.5 kg/m^3','0.7 kg/m^3','1 kg/m^3',...
       '2 kg/m^3','3 kg/m^3','5 kg/m^3','7 kg/m^3','10 kg/m^3','20 kg/m^3','50 kg/m^3','100 kg/m^3' ,...
       'Location','BestOutside')
print -dpng Multiple_cases_noset_C_D.png
%%
figure(125);clf
ylabel('settling velocity (mm s^{-1})','FontSize',14)
xlabel('Average Concentration (kg m^{-3})','FontSize',14)
a14=line(c{14}, wsa{14},'Marker','.','MarkerSize',16,'Color',[.5,0,.5],'LineStyle','none');%100 %77
a13=line(c{13}, wsa{13},'Marker','.','MarkerSize',16,'Color',[0,0,.5],'LineStyle','none');%50   %76
a12=line(c{12}, wsa{12},'Marker','.','MarkerSize',16,'Color',[0,.5,.5],'LineStyle','none');%20  %75
a11=line(c{11}, wsa{11},'Marker','.','MarkerSize',16,'Color',[.8,.8,.8],'LineStyle','none');%10 %74
a10=line(c{10}, wsa{10},'Marker','.','MarkerSize',16,'Color','y','LineStyle','none');%7         %73
a3=line(c{3}, wsa{3},'Marker','.','MarkerSize',16,'Color',[.5,0,0],'LineStyle','none');%5       %58
a9=line(c{9}, wsa{9},'Marker','.','MarkerSize',16,'Color','r','LineStyle','none');%3            %72
a8=line(c{8}, wsa{8},'Marker','.','MarkerSize',16,'Color','c','LineStyle','none');%2            %71
a1=line(c{1}, wsa{1},'Marker','.','MarkerSize',16,'Color',[.8,.5,0],'LineStyle','none');%1      %56
a7=line(c{7}, wsa{7},'Marker','.','MarkerSize',16,'Color','k','LineStyle','none');%0.7          %70
a6=line(c{6}, wsa{6},'Marker','.','MarkerSize',16,'Color','b','LineStyle','none');%0.5          %69
a5=line(c{5}, wsa{5},'Marker','.','MarkerSize',16,'Color','g','LineStyle','none');%0.3          %68
a2=line(c{2}, wsa{2},'Marker','.','MarkerSize',16,'Color',[.3,.3,.3],'LineStyle','none');%0.2   %57
a4=line(c{4}, wsa{4},'Marker','.','MarkerSize',16,'Color',[.6,.6,.6],'LineStyle','none');%0.1   %67
axis([0, 10,0,6])
legend([a4(1),a2(1),a5(1),a6(1),a7(1),a1(1),a8(1),a9(1),a3(1),a10(1),a11(1),a12(1),a13(1),a14(1)],...
       '0.1 kg/m^3','0.2 kg/m^3','0.3 kg/m^3','0.5 kg/m^3','0.7 kg/m^3','1 kg/m^3',...
       '2 kg/m^3','3 kg/m^3','5 kg/m^3','7 kg/m^3','10 kg/m^3','20 kg/m^3','50 kg/m^3','100 kg/m^3' ,...
       'Location','BestOutside')
print -dpng Multiple_cases_noset_C_WS.png
%%
figure(126);clf
ylabel('settling velocity (mm s^{-1})','FontSize',14)
xlabel('Diameter D (\mum)','FontSize',14)
a14=line(d{14}, wsa{14},'Marker','.','MarkerSize',16,'Color',[.5,0,.5],'LineStyle','none');%100 %77
a13=line(d{13}, wsa{13},'Marker','.','MarkerSize',16,'Color',[0,0,.5],'LineStyle','none');%50   %76
a12=line(d{12}, wsa{12},'Marker','.','MarkerSize',16,'Color',[0,.5,.5],'LineStyle','none');%20  %75
a11=line(d{11}, wsa{11},'Marker','.','MarkerSize',16,'Color',[.8,.8,.8],'LineStyle','none');%10 %74
a10=line(d{10}, wsa{10},'Marker','.','MarkerSize',16,'Color','y','LineStyle','none');%7         %73
a3=line(d{3}, wsa{3},'Marker','.','MarkerSize',16,'Color',[.5,0,0],'LineStyle','none');%5       %58
a9=line(d{9}, wsa{9},'Marker','.','MarkerSize',16,'Color','r','LineStyle','none');%3            %72
a8=line(d{8}, wsa{8},'Marker','.','MarkerSize',16,'Color','c','LineStyle','none');%2            %71
a1=line(d{1}, wsa{1},'Marker','.','MarkerSize',16,'Color',[.8,.5,0],'LineStyle','none');%1      %56
a7=line(d{7}, wsa{7},'Marker','.','MarkerSize',16,'Color','k','LineStyle','none');%0.7          %70
a6=line(d{6}, wsa{6},'Marker','.','MarkerSize',16,'Color','b','LineStyle','none');%0.5          %69
a5=line(d{5}, wsa{5},'Marker','.','MarkerSize',16,'Color','g','LineStyle','none');%0.3          %68
a2=line(d{2}, wsa{2},'Marker','.','MarkerSize',16,'Color',[.3,.3,.3],'LineStyle','none');%0.2   %57
a4=line(d{4}, wsa{4},'Marker','.','MarkerSize',16,'Color',[.6,.6,.6],'LineStyle','none');%0.1   %67
%axis([0, 10,0,6])
legend([a4(1),a2(1),a5(1),a6(1),a7(1),a1(1),a8(1),a9(1),a3(1),a10(1),a11(1),a12(1),a13(1),a14(1)],...
       '0.1 kg/m^3','0.2 kg/m^3','0.3 kg/m^3','0.5 kg/m^3','0.7 kg/m^3','1 kg/m^3',...
       '2 kg/m^3','3 kg/m^3','5 kg/m^3','7 kg/m^3','10 kg/m^3','20 kg/m^3','50 kg/m^3','100 kg/m^3' ,...
       'Location','BestOutside')
print -dpng Multiple_cases_noset_D_WS.png
%%
figure(127);clf
ylabel('settling velocity (mm s^{-1})','FontSize',14)
xlabel('Turbulence Shear Rate G (s^{-1})','FontSize',14)
a14=line(G{14}, wsa{14},'Marker','.','MarkerSize',16,'Color',[.5,0,.5],'LineStyle','none');%100 %77
a13=line(G{13}, wsa{13},'Marker','.','MarkerSize',16,'Color',[0,0,.5],'LineStyle','none');%50   %76
a12=line(G{12}, wsa{12},'Marker','.','MarkerSize',16,'Color',[0,.5,.5],'LineStyle','none');%20  %75
a11=line(G{11}, wsa{11},'Marker','.','MarkerSize',16,'Color',[.8,.8,.8],'LineStyle','none');%10 %74
a10=line(G{10}, wsa{10},'Marker','.','MarkerSize',16,'Color','y','LineStyle','none');%7         %73
a3=line(G{3}, wsa{3},'Marker','.','MarkerSize',16,'Color',[.5,0,0],'LineStyle','none');%5       %58
a9=line(G{9}, wsa{9},'Marker','.','MarkerSize',16,'Color','r','LineStyle','none');%3            %72
a8=line(G{8}, wsa{8},'Marker','.','MarkerSize',16,'Color','c','LineStyle','none');%2            %71
a1=line(G{1}, wsa{1},'Marker','.','MarkerSize',16,'Color',[.8,.5,0],'LineStyle','none');%1      %56
a7=line(G{7}, wsa{7},'Marker','.','MarkerSize',16,'Color','k','LineStyle','none');%0.7          %70
a6=line(G{6}, wsa{6},'Marker','.','MarkerSize',16,'Color','b','LineStyle','none');%0.5          %69
a5=line(G{5}, wsa{5},'Marker','.','MarkerSize',16,'Color','g','LineStyle','none');%0.3          %68
a2=line(G{2}, wsa{2},'Marker','.','MarkerSize',16,'Color',[.3,.3,.3],'LineStyle','none');%0.2   %57
a4=line(G{4}, wsa{4},'Marker','.','MarkerSize',16,'Color',[.6,.6,.6],'LineStyle','none');%0.1   %67
axis([0, 10,0,6])
legend([a4(1),a2(1),a5(1),a6(1),a7(1),a1(1),a8(1),a9(1),a3(1),a10(1),a11(1),a12(1),a13(1),a14(1)],...
       '0.1 kg/m^3','0.2 kg/m^3','0.3 kg/m^3','0.5 kg/m^3','0.7 kg/m^3','1 kg/m^3',...
       '2 kg/m^3','3 kg/m^3','5 kg/m^3','7 kg/m^3','10 kg/m^3','20 kg/m^3','50 kg/m^3','100 kg/m^3' ,...
       'Location','BestOutside')
print -dpng Multiple_cases_noset_G_WS.png



%% plot all results
figure(122);clf
ylabel('Average Diameter D (\mum)','FontSize',14)
xlabel('Turbulence Shear Rate G (s^{-1})','FontSize',14)
a14=line(Gall{14}, dall{14},'Marker','.','MarkerSize',16,'Color',[.5,0,.5],'LineStyle','none');%100 %77
a13=line(Gall{13}, dall{13},'Marker','.','MarkerSize',16,'Color',[0,0,.5],'LineStyle','none');%50   %76
a12=line(Gall{12}, dall{12},'Marker','.','MarkerSize',16,'Color',[0,.5,.5],'LineStyle','none');%20  %75
a11=line(Gall{11}, dall{11},'Marker','.','MarkerSize',16,'Color',[.8,.8,.8],'LineStyle','none');%10 %74
a10=line(Gall{10}, dall{10},'Marker','.','MarkerSize',16,'Color','y','LineStyle','none');%7         %73
a3=line(Gall{3}, dall{3},'Marker','.','MarkerSize',16,'Color',[.5,0,0],'LineStyle','none');%5       %58
a9=line(Gall{9}, dall{9},'Marker','.','MarkerSize',16,'Color','r','LineStyle','none');%3            %72
a8=line(Gall{8}, dall{8},'Marker','.','MarkerSize',16,'Color','c','LineStyle','none');%2            %71
a1=line(Gall{1}, dall{1},'Marker','.','MarkerSize',16,'Color',[.8,.5,0],'LineStyle','none');%1      %56
a7=line(Gall{7}, dall{7},'Marker','.','MarkerSize',16,'Color','k','LineStyle','none');%0.7          %70
a6=line(Gall{6}, dall{6},'Marker','.','MarkerSize',16,'Color','b','LineStyle','none');%0.5          %69
a5=line(Gall{5}, dall{5},'Marker','.','MarkerSize',16,'Color','g','LineStyle','none');%0.3          %68
a2=line(Gall{2}, dall{2},'Marker','.','MarkerSize',16,'Color',[.3,.3,.3],'LineStyle','none');%0.2   %57
a4=line(Gall{4}, dall{4},'Marker','.','MarkerSize',16,'Color',[.6,.6,.6],'LineStyle','none');%0.1   %67
%axis([0, 6,100,700])
legend([a4(1),a2(1),a5(1),a6(1),a7(1),a1(1),a8(1),a9(1),a3(1),a10(1),a11(1),a12(1),a13(1),a14(1)],...
       '0.1 kg/m^3','0.2 kg/m^3','0.3 kg/m^3','0.5 kg/m^3','0.7 kg/m^3','1 kg/m^3',...
       '2 kg/m^3','3 kg/m^3','5 kg/m^3','7 kg/m^3','10 kg/m^3','20 kg/m^3','50 kg/m^3','100 kg/m^3' ,...
       'Location','BestOutside')
print -dpng Multiple_cases_all_noset.png


%%
return
%%

nf = 2.;
m = .2;
q = (nf-1.)/(2.*m);
q19 = (1.9-1.)/(2.*m);
q21 = (2.1-1.)/(2.*m);

figure(13);clf
%plot(tconc_mn/Gc_mn**q,(1e-3*fdiam_av_mn)**(2*q),'.')
%plot( c[0]/G[0]**q, (1e-3*d[0])**(2.*q),'.',label="1 kg/m^3")
line(c{4}./G{4}.^q, d{4}.^(2*q),'Marker','.','MarkerSize',16,'Color',[.6,.6,.6],'LineStyle','none')%0.1   %67
line(c{2}./G{2}.^q, d{2}.^(2*q),'Marker','.','MarkerSize',16,'Color',[.3,.3,.3],'LineStyle','none')%0.2   %57
line(c{5}./G{5}.^q, d{5}.^(2*q),'Marker','.','MarkerSize',16,'Color','g','LineStyle','none')%0.3          %68
line(c{6}./G{6}.^q, d{6}.^(2*q),'Marker','.','MarkerSize',16,'Color','b','LineStyle','none')%0.5          %69
line(c{7}./G{7}.^q, d{7}.^(2*q),'Marker','.','MarkerSize',16,'Color','k','LineStyle','none')%0.7          %70
line(c{1}./G{1}.^q, d{1}.^(2*q),'Marker','.','MarkerSize',16,'Color',[.8,.5,0],'LineStyle','none')%1      %56
line(c{8}./G{8}.^q, d{8}.^(2*q),'Marker','.','MarkerSize',16,'Color','c','LineStyle','none')%2            %71
line(c{9}./G{9}.^q, d{9}.^(2*q),'Marker','.','MarkerSize',16,'Color','r','LineStyle','none')%3            %72
line(c{3}./G{3}.^q, d{3}.^(2*q),'Marker','.','MarkerSize',16,'Color',[.5,0,0],'LineStyle','none')%5       %58
line(c{10}./G{10}.^q, d{10}.^(2*q),'Marker','.','MarkerSize',16,'Color','y','LineStyle','none')%7         %73
line(c{11}./G{11}.^q, d{11}.^(2*q),'Marker','.','MarkerSize',16,'Color',[.8,.8,.8],'LineStyle','none')%10 %74
line(c{12}./G{12}.^q, d{12}.^(2*q),'Marker','.','MarkerSize',16,'Color',[0,.5,.5],'LineStyle','none')%20  %75
line(c{13}./G{13}.^q, d{13}.^(2*q),'Marker','.','MarkerSize',16,'Color',[0,0,.5],'LineStyle','none')%50   %76
line(c{14}./G{14}.^q, d{14}.^(2*q),'Marker','.','MarkerSize',16,'Color',[.5,0,.5],'LineStyle','none')%100 %77
title(['Winterwerp plot: nf=',num2str(nf),' m=',num2str(m),' q=',num2str(q)])
ylabel('Diameter ^(2q)')
xlabel('Conc/G^q')
legend('0.1 kg/m^3','0.2 kg/m^3','0.3 kg/m^3','0.5 kg/m^3','0.7 kg/m^3','1 kg/m^3',...
       '2 kg/m^3','3 kg/m^3','5 kg/m^3','7 kg/m^3','10 kg/m^3','20 kg/m^3','50 kg/m^3','100 kg/m^3' ,...
       'Location','BestOutside')
print -dpng Winterwerp_Multiple_cases_noset.png
%%
figure(133);clf
%plot(tconc_mn/Gc_mn**q,(1e-3*fdiam_av_mn)**(2*q),'.')
%plot( c[0]/G[0]**q, (1e-3*d[0])**(2.*q),'.',label="1 kg/m^3")
a14=line(call{14}./Gall{14}.^q, dall{14}.^(2*q),'Marker','.','MarkerSize',16,'Color',[.5,0,.5],'LineStyle','none');%100 %77
a13=line(call{13}./Gall{13}.^q, dall{13}.^(2*q),'Marker','.','MarkerSize',16,'Color',[0,0,.5],'LineStyle','none');%50   %76
a12=line(call{12}./Gall{12}.^q, dall{12}.^(2*q),'Marker','.','MarkerSize',16,'Color',[0,.5,.5],'LineStyle','none');%20  %75
a11=line(call{11}./Gall{11}.^q, dall{11}.^(2*q),'Marker','.','MarkerSize',16,'Color',[.8,.8,.8],'LineStyle','none');%10 %74
a10=line(call{10}./Gall{10}.^q, dall{10}.^(2*q),'Marker','.','MarkerSize',16,'Color','y','LineStyle','none');%7         %73
a3=line(call{3}./Gall{3}.^q, dall{3}.^(2*q),'Marker','.','MarkerSize',16,'Color',[.5,0,0],'LineStyle','none');%5       %58
a9=line(call{9}./Gall{9}.^q, dall{9}.^(2*q),'Marker','.','MarkerSize',16,'Color','r','LineStyle','none');%3            %72
a8=line(call{8}./Gall{8}.^q, dall{8}.^(2*q),'Marker','.','MarkerSize',16,'Color','c','LineStyle','none');%2            %71
a1=line(call{1}./Gall{1}.^q, dall{1}.^(2*q),'Marker','.','MarkerSize',16,'Color',[.8,.5,0],'LineStyle','none');%1      %56
a7=line(call{7}./Gall{7}.^q, dall{7}.^(2*q),'Marker','.','MarkerSize',16,'Color','k','LineStyle','none');%0.7          %70
a6=line(call{6}./Gall{6}.^q, dall{6}.^(2*q),'Marker','.','MarkerSize',16,'Color','b','LineStyle','none');%0.5          %69
a5=line(call{5}./Gall{5}.^q, dall{5}.^(2*q),'Marker','.','MarkerSize',16,'Color','g','LineStyle','none');%0.3          %68
a2=line(call{2}./Gall{2}.^q, dall{2}.^(2*q),'Marker','.','MarkerSize',16,'Color',[.3,.3,.3],'LineStyle','none');%0.2   %57
a4=line(call{4}./Gall{4}.^q, dall{4}.^(2*q),'Marker','.','MarkerSize',16,'Color',[.6,.6,.6],'LineStyle','none');%0.1   %67
title(['Winterwerp plot: nf=',num2str(nf),' m=',num2str(m),' q=',num2str(q)])
ylabel('Diameter ^(2q)')
xlabel('Conc/G^q')
legend([a4(1),a2(1),a5(1),a6(1),a7(1),a1(1),a8(1),a9(1),a3(1),a10(1),a11(1),a12(1),a13(1),a14(1)],...
       '0.1 kg/m^3','0.2 kg/m^3','0.3 kg/m^3','0.5 kg/m^3','0.7 kg/m^3','1 kg/m^3',...
       '2 kg/m^3','3 kg/m^3','5 kg/m^3','7 kg/m^3','10 kg/m^3','20 kg/m^3','50 kg/m^3','100 kg/m^3' ,...
       'Location','BestOutside')
print -dpng Winterwerp_Multiple_cases_all_noset.png
%%
figure(14);clf
xlabel('Average concentration ','FontSize',14)
ylabel('settling velocity (s^{-1})','FontSize',14)
line(c{4},wsa{4},'Marker','.','MarkerSize',16,'Color',[.6,.6,.6],'LineStyle','none')%0.1   %67
line(c{2},wsa{2},'Marker','.','MarkerSize',16,'Color',[.3,.3,.3],'LineStyle','none')%0.2   %57
line(c{5},wsa{5},'Marker','.','MarkerSize',16,'Color','g','LineStyle','none')%0.3          %68
line(c{6},wsa{6},'Marker','.','MarkerSize',16,'Color','b','LineStyle','none')%0.5          %69
line(c{7},wsa{7},'Marker','.','MarkerSize',16,'Color','k','LineStyle','none')%0.7          %70
line(c{1},wsa{1},'Marker','.','MarkerSize',16,'Color',[.8,.5,0],'LineStyle','none')%1      %56
line(c{8},wsa{8},'Marker','.','MarkerSize',16,'Color','c','LineStyle','none')%2            %71
line(c{9},wsa{9},'Marker','.','MarkerSize',16,'Color','r','LineStyle','none')%3            %72
line(c{3},wsa{3},'Marker','.','MarkerSize',16,'Color',[.5,0,0],'LineStyle','none')%5       %58
line(c{10},wsa{10},'Marker','.','MarkerSize',16,'Color','y','LineStyle','none')%7         %73
line(c{11},wsa{11},'Marker','.','MarkerSize',16,'Color',[.8,.8,.8],'LineStyle','none')%10 %74
line(c{12},wsa{12},'Marker','.','MarkerSize',16,'Color',[0,.5,.5],'LineStyle','none')%20  %75
line(c{13},wsa{13},'Marker','.','MarkerSize',16,'Color',[0,0,.5],'LineStyle','none')%50   %76
line(c{14},wsa{14},'Marker','.','MarkerSize',16,'Color',[.5,0,.5],'LineStyle','none')%100 %77
set(gca,'YScale','log','XScale','log');grid
%axis([0, 6,100,700])
legend('0.1 kg/m^3','0.2 kg/m^3','0.3 kg/m^3','0.5 kg/m^3','0.7 kg/m^3','1 kg/m^3',...
       '2 kg/m^3','3 kg/m^3','5 kg/m^3','7 kg/m^3','10 kg/m^3','20 kg/m^3','50 kg/m^3','100 kg/m^3' ,...
       'Location','BestOutside')
print -dpng Multiple_cases_diam_WS_noset.png

