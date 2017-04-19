clear

% Winterwerp
% wD50=wDp+ka*C0/(kb*sqrt(G))

flist = {'mor01','mor02','mor03','mor05','mor07','mor1',...
         'mor2','mor3','mor5','mor7','mor10','mor20','mor50','mor100'};
for ii=1:length(flist)
    run_name=flist{ii};
    G{ii} =load([run_name,'gc.txt']);
    d{ii}=load([run_name,'diam.txt']);
    wsa{ii}=load([run_name,'ws.txt']);
    c{ii}=load([run_name,'mconc.txt']);
end
%%
figure(12);clf
ylabel('Average Diameter D (\mum)','FontSize',14)
xlabel('C/sqrt(G) ','FontSize',14)
line(c{1}./sqrt(G{1}),d{1},'Marker','.','Color',[.6,.6,.6],'LineStyle','none')
line(c{2}./sqrt(G{2}),d{2},'Marker','.','Color',[.3,.3,.3],'LineStyle','none')
line(c{3}./sqrt(G{3}),d{3},'Marker','.','Color','g','LineStyle','none')
line(c{4}./sqrt(G{4}),d{4},'Marker','.','Color','b','LineStyle','none')
line(c{5}./sqrt(G{5}), d{5},'Marker','.','Color','k','LineStyle','none')
line(c{6}./sqrt(G{6}), d{6},'Marker','.','Color',[.8,.5,0],'LineStyle','none')
line(c{7}./sqrt(G{7}), d{7},'Marker','.','Color','c','LineStyle','none')
line(c{8}./sqrt(G{8}), d{8},'Marker','.','Color','r','LineStyle','none')
line(c{9}./sqrt(G{9}), d{9},'Marker','.','Color',[.5,0,0],'LineStyle','none')
line(c{10}./sqrt(G{10}), d{10},'Marker','.','Color','y','LineStyle','none')
line(c{11}./sqrt(G{11}), d{11},'Marker','.','Color',[.8,.8,.8],'LineStyle','none')
line(c{12}./sqrt(G{12}), d{12},'Marker','.','Color',[0,.5,.5],'LineStyle','none')
line(c{13}./sqrt(G{13}), d{13},'Marker','.','Color',[0,0,.5],'LineStyle','none')
%line(c{14}./sqrt(G{14}), d{14},'Marker','.','Color',[.5,0,.5],'LineStyle','none')
axis([0, 2,0,3500])
legend('0.1 kg/m^3','0.2 kg/m^3','0.3 kg/m^3','0.5 kg/m^3','0.7 kg/m^3',...
       '1 kg/m^3','2 kg/m^3' ,'3 kg/m^3','5 kg/m^3','7 kg/m^3',...
       '10 kg/m^3','20 kg/m^3','50 kg/m^3','100 kg/m^3','Location','BestOutside')
%%
dal=[];
cgal=[];
for ii=1:length(flist)
    ia=find(d{ii}<3000&c{ii}./sqrt(G{ii})<1);
    dal=[dal,d{ii}(ia)];
    cgal=[cgal,c{ii}(ia)./sqrt(G{ii}(ia))];
end
[cgal,I]=sort(cgal);
dal=dal(I)/1e6;
%%
[taub,tau,h,sig,Z,S,sigma,sen,n,senplot,CIlower,CIupper,D,Dall,C3]= ktaub([cgal',dal'], .05, 0);
myFit = LinearModel.fit(cgal,dal);inter=myFit.Coefficients.Estimate(1);
B=robustfit(cgal,dal);
X1 = ones(length(cgal),2); 
X1(:,2) = cgal';
B1 = X1\dal'; % 

figure(13);clf
ylabel('Average Diameter D (\mum)','FontSize',14)
xlabel('C/sqrt(G) ','FontSize',14)
line(cgal,dal,'Marker','.','Color',[.1,.1,.1],'LineStyle','none')
axis([0, 1,0,3000/1e6])
hold on
line([0:.1:1],B(1)+[0:.1:1]*B(2),'Color','r')
h0=shadedErrorBar([0:.1:1],[0:.1:1]*sen+inter,[[0:.1:1]'*(CIupper-sen),[0:.1:1]'*(-CIlower+sen)]','b');
text(.1,2.5e-3,['D_{50}=',num2str(sen*1e3,'%3.2f'),' 10^{-3} C/sqrt(G)+',num2str(inter*1e5,'%3.1f'),' 10^{-5}'])
text(.2,2e-3,['K_a/K_b=',num2str(sen*1e3,'%3.2f'),' 10^{-3}'],'FontWeight','bold')

print -dpng -painters winterwerp_model_comp.png
%%

for ii=1:length(flist)-1
    run_name=flist{ii};
    load([run_name,'_all.mat'],'Gc');
    Gall{ii} =Gc(10:end,:);
    load([run_name,'_all.mat'],'fdiam_av');
    dall{ii}=fdiam_av(10:end,:);
    load([run_name,'_all.mat'],'ws_av');
    wsaall{ii}=ws_av(10:end,:);
    load([run_name,'_all.mat'],'tconc');
    call{ii}=tconc(10:end,:);
end

dal=[];
cgal=[];
for ii=1:length(flist)-1
    ia=find(dall{ii}<3000&call{ii}./sqrt(Gall{ii})<1);
    dal=[dal,dall{ii}(ia)'];
    cgal=[cgal,call{ii}(ia)'./sqrt(Gall{ii}(ia)')];
end
[cgal,I]=sort(cgal);
dal=dal(I)/1e6;
%%
%[taub,tau,h,sig,Z,S,sigma,sen,n,senplot,CIlower,CIupper,D,Dall,C3]= ktaub([cgal',dal'], .05, 0);
myFit = LinearModel.fit(cgal,dal);inter=myFit.Coefficients.Estimate(1);
B=robustfit(cgal,dal);
X1 = ones(length(cgal),2); 
X1(:,2) = cgal';
B1 = X1\dal'; % 
sen=B1(2);

figure(14);clf
ylabel('Average Diameter D (\mum)','FontSize',14)
xlabel('C/sqrt(G) ','FontSize',14)
line(cgal,dal,'Marker','.','Color',[.1,.1,.1],'LineStyle','none')
axis([0, 1,0,3000/1e6])
hold on
line([0:.1:1],B(1)+[0:.1:1]*B(2),'Color','r')
h0=shadedErrorBar([0:.1:1],[0:.1:1]*sen+inter,[[0:.1:1]'*(CIupper-sen),[0:.1:1]'*(-CIlower+sen)]','b');
text(.1,2.5e-3,['D_{50}=',num2str(sen*1e3,'%3.2f'),' 10^{-3} C/sqrt(G)+',num2str(inter*1e5,'%3.1f'),' 10^{-5}'])
text(.2,2e-3,['K_a/K_b=',num2str(sen*1e3,'%3.2f'),' 10^{-3}'],'FontWeight','bold')

