clear

I = imread('Verneyetal_CSR_2011_Fig7_crop.jpg');

figure(1);clf
set(gcf,'PaperPosition',[.5,.5,10,6]);wysiwyg
image(I);
set(gca,'Visible','off')
aa0=get(gca,'Position');aa=aa0;
aa(1)=0.04;aa(3)=0.8;
set(gca,'Position',aa)
%%
load Verney_multiinit.mat 
url='http://geoport.whoi.edu/thredds/dodsC/sand/usgs/users/aretxabaleta/floc/ocean_his_';
nc = ncgeodataset([url,'verneyc13.nc'])
fdiam = 1e6*nc{'Sd50'}(:)

%%
lineStyles=linspecer(length(allr));

%%
axes('Position',[0.0860    0.1750    0.7410    0.7410])
set(gca,'Visible','off')
hold on
for ii=1:length(allr)
    h1(ii)=line(allr{ii}.tt,allr{ii}.fdiam_25,'Color',lineStyles(ii,:),'LineWidth',3.5);
end
%aa=get(gca,'Position');
set(gca,'YLim',[40,1500],'XLim',[0,700],'YScale','log','FontSize',14,'FontWeight','bold')
ylabel('Average diameter (mm) ')
xlabel('Time (min)');grid
hle=legend([h1],['C_0',num2str(fdiam(1),'%2.0f'),'\mum'],['C_0',num2str(fdiam(3),'%2.0f'),'\mum'],...
    ['C_0',num2str(fdiam(5),'%2.0f'),'\mum'],['C_0',num2str(fdiam(7),'%2.0f'),'\mum'],...
    ['C_0',num2str(fdiam(9),'%2.0f'),'\mum'],['C_0',num2str(fdiam(11),'%2.0f'),'\mum'],...
    ['C_0',num2str(fdiam(13),'%2.0f'),'\mum'],['C_0',num2str(fdiam(15),'%2.0f'),'\mum'],...
    ['C_0',num2str(fdiam(17),'%2.0f'),'\mum'],['C_0',num2str(fdiam(19),'%2.0f'),'\mum'],...
    ['C_0',num2str(fdiam(21),'%2.0f'),'\mum'],...
       'Location','BestOutside');
set(gca,'Position',[0.0860    0.1750    0.7410    0.7410])

print -dpng -painters compar_figure7_verney.png