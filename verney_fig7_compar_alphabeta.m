clear

I = imread('Verneyetal_CSR_2011_Fig7_crop.jpg');

load('verney_obs.mat');

figure(1);clf
set(gcf,'PaperPosition',[.5,.5,10,6]);wysiwyg
image(I);
set(gca,'Visible','off')
aa0=get(gca,'Position');aa=aa0;
aa(1)=0.03;aa(3)=0.8;
set(gca,'Position',aa)
%%
load Verney_multiinit_alphabeta.mat 
url='http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/aretxabaleta/floc/ocean_his_';
nc = ncgeodataset([url,'verneyc11_ncs15.nc'])
fdiam = 1e6*nc{'Sd50'}(:);

%%
lineStyles=linspecer(length(allr));

%%
axes('Position',[0.0760    0.1750    0.7410    0.7410])
set(gca,'Visible','off')
hold on
for ii=1:length(allr)
    %h1(ii)=line(allr{ii}.tt,allr{ii}.fdiam_25,'Color',lineStyles(ii,:),'LineWidth',3.5);
    modi=interp1(allr{ii}.tt,allr{ii}.fdiam_25,ttobso);
    err(ii)=rms(modi-diamobso);    
end
for ii=1:length(allr)
    %h1(ii)=line(allr{ii}.tt,allr{ii}.fdiam_ar,'Color',lineStyles(ii,:),'LineWidth',3.5);
    modia=interp1(allr{ii}.tt,allr{ii}.fdiam_ar,ttobso);
    erra(ii)=rms(modia-diamobso);    
end
for ii=1:length(allr)
    h1(ii)=line(allr{ii}.tt,allr{ii}.fdiam_dum,'Color',lineStyles(ii,:),'LineWidth',3.5);
    modia3=interp1(allr{ii}.tt,allr{ii}.fdiam_dum,ttobso);
    erra3(ii)=rms(modia3-diamobso);    
end
%aa=get(gca,'Position');
%line(ttobs,diamobs,'LineStyle','none','Marker','.','MarkerFaceColor','b','MarkerSize',20)
for ii=1:length(allr)
    text(allr{ii}.tt(37),allr{ii}.fdiam_dum(37),...
        ['rmse=',num2str(erra3(ii),'%2.0f'),'\mum']);
end
set(gca,'YLim',[40,1500],'XLim',[0,700],'YScale','log','FontSize',14,'FontWeight','bold')
ylabel('Average diameter (mm) ')
xlabel('Time (min)');grid
hle=legend([h1],'\alpha 0.35 \beta 0.12','\alpha 0.25 \beta 0.1','\alpha 0.35 \beta 0.08',...
    '\alpha 0.40 \beta 0.12','\alpha 0.25 \beta 0.12','\alpha 0.35 \beta 0.1',...
    '\alpha 0.50 \beta 0.12','\alpha 0.25 \beta 0.15','\alpha 0.35 \beta 0.14','\alpha 0.50 \beta 0.15',...
    '\alpha 0.30 \beta 0.12','\alpha 0.35 \beta 0.20',...
       'Location','BestOutside');
set(gca,'Position',[0.0760    0.1750    0.7410    0.7410])

print -dpng -painters compar_figure7_verney_alphabeta.png

%%
row_text={[],[],[],[],[],[],[],[],[],[],[],[]};
col_text={'\alpha','\beta','RMSE'};
D=[0.35,0.25,0.35,0.4,0.25,0.35,0.5,0.25,0.35,0.5,0.3,0.35;...
   0.12,0.1,0.08,0.12,0.12,0.1,0.12,0.15,0.14,0.15,0.12,0.2;...
   erra3]';
errtab=html_table(D,row_text,col_text,['RMSE'],'alphabeta_rms.html',1)

