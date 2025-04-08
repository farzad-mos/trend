load('F:\Analysis\BS\trend\data\sa\all\SA_1.mat')


%% remove gross errors & outliers 

clearvars  h2 a1 a2 e1 e2 h1


%convert SSH to DT
sa.dt=sa.ssh-sa.N+sa.dac+sa.Ell_corr;


%   find gross errors
sa.ssh(abs(sa.dt)>=2)=NaN;

%   find random errors
for i=min(sa.said):max(sa.said)
    
    id=unique(sa.id(sa.said==i));
    if ~isempty(id)
        for j=id(1):id(end)
            
            if ~isempty(sa.ssh(sa.said==i&sa.id==j))
                
                m=mean(abs(sa.dt(sa.said==i&sa.id==j)),'omitnan');
                s=std(abs(sa.dt(sa.said==i&sa.id==j)),'omitnan');
                sa.ssh(abs(sa.dt(sa.said==i&sa.id==j))>=m+(s*3))=NaN;
                sa.ssh(abs(sa.dt(sa.said==i&sa.id==j))<=m-(s*3))=NaN;

            end
        end
    end
end

sa = rmmissing(sa,'DataVariables',{'ssh'}); %remove NaN data

% find outliers
for k=1:5

%   remove outlier data using scaled MAD
sacor=table();

for i=min(sa.said):max(sa.said)
    id=unique(sa.id(sa.said==i));
    if ~isempty(id)
        for j=id(1):id(end)
            if ~isempty(sa.ssh(sa.said==i&sa.id==j))
                sacor1=sa(sa.said==i&sa.id==j,:);
                [~,out]=rmoutliers(sacor1.dt,'movmedian',100);
                [m, ~]=find(out==1);
                sacor1.ssh(m)=NaN;
                sacor=[sacor;sacor1];
                clear sacor1
            end
        end
    end
end

sacor = rmmissing(sacor,'DataVariables',{'ssh'}); %remove NaN data
clear sa
sa=sacor;
clear sacor

end

sa.year_sa=year(sa.t);

clearvars i k 


%% find within 50km distance data at TG
% extract points within 50 km distance of TG

for i=1:height(TGinfo)
[lat,lon] = scircle1(TGinfo.Lat(i),TGinfo.Lon(i),km2deg(50),[],[],[],20);
in = inpolygon(sa.lat,sa.lon,lat,lon);
sa_tg{i}=sa(in==1,:);
% [~,d] = dsearchn([Lat Lon],[sa.lat sa.lon]);
clearvars i in lat lon
end
clear i

%% TG treatment
% combine all TG tables
tg={TG2,TG3,TG4,TG6,TG7,TG8,TG9,TG10,TG11,TG12,TG13,TG14};
% load('F:\Analysis\BS\trend\ASL trend\TG_4st.mat') %LU added 

% tide system
NAPlat=dms2degrees([52 22 53]);
% for i=1:height(TGinfo)


% T=eval(['TG' num2str(i)]);
% [~,ia] = unique(T.date);
% tg{i}=T(ia,:);


for i=1:length(tg)
%TIDE system correction; zero to mean
tg{i}.zerotomean=repmat(0.29541*(sind(TGinfo.Lat(i)).^2-sind(NAPlat).^2)+0.00042*(sind(TGinfo.Lat(i)).^4-sind(NAPlat).^4),height(tg{i}),1).*100; %add to TG data
tg{i}.dt2=tg{i}.dt+tg{i}.zerotomean;

clearvars T ia
end
clearvars i NAPlat

%% SA treatment

% combine all sa
% load sa_cut_all

for i=1:15
sa_tg{i}=[sa1{i};sa2{i};sa3{i};sa4{i};sa5{i};sa6{i};sa7{i};sa8{i};sa9{i}];
sa_tg{i} = sortrows(sa_tg{i},'t','ascend');
end

% Ellipsoid correction

% t/p
a1=6378136.3;
e1=0.081819221456;
% GRS80

a2=6378137.0;
e2=0.081819190842621;

for i=1:15
h1=(a1*(1-e1^2))./(1-(e1^2)*((sind(sa_tg{i}.lat)).^2)).^.5;
h2=(a2*(1-e2^2))./(1-(e2^2)*((sind(sa_tg{i}.lat)).^2)).^.5;
Ell_corr=h1-h2;

sa_tg{i}.dt = (sa_tg{i}.dt+Ell_corr);
clear Ell_corr
end
clearvars a1 a2 e1 e2 h1 h2 sa1 sa2 sa3 sa4 sa5 sa6 sa7 sa8 sa9 Ell_corr i 

%% 

for i=1:15
sa_tg{i} = removevars(sa_tg{i}, {'dac','oct'});


[m1,~]=find(sa10{i}.t>=max(sa8{i}.t));
[m2,~]=find(sa11{i}.t>=max(sa9{i}.t));


sa{i}=[sa_tg{i};sa10{i}(m1,:);sa11{i}(m2,:)];
sa{i} = sortrows(sa{i},'t','ascend');

clearvars m1 m2
end

%% mean SA at eact TG per cycle
%  load('F:\Analysis\BS\trend\ASL trend\sa_tg.mat')


for k=1:height(TGinfo)
    sacor=table();
    H=1;
    for i=min(sa_tg{k}.said):max(sa_tg{k}.said)
        
        id=unique(sa_tg{1, k}.id(sa_tg{1, k}.said==i));
        if ~isempty(id)
            for j=id(1):id(end)
                
                if ~isempty(sa_tg{1, k}.ssh(sa_tg{1, k}.said==i&sa_tg{1, k}.id==j))
                    
                    sacor.dt(H)=mean((sa_tg{1, k}.dt(sa_tg{1, k}.said==i&sa_tg{1, k}.id==j)),'omitnan');
                    sacor.dac(H)=mean((sa_tg{1, k}.dac(sa_tg{1, k}.said==i&sa_tg{1, k}.id==j)),'omitnan');
                    time=sa_tg{1, k}.t(sa_tg{1, k}.said==i&sa_tg{1, k}.id==j);
                    sacor.t(H)=time(1);
                    sacor.said(H)=i;
%                   sacor.tgdt(H)=interp1(tg{1, k}.Time,tg{1, k}.dt,time(1));
                    H=H+1;
                    clear time
                end
            end
        end
    end
    satg{k}=sacor;
    clear sacor
end

%% data preparation (monthly mean)
sacut=satg;

% sort by date
for k=1:height(TGinfo)
sacut{k} = sortrows(sacut{k},'t','ascend');
% add monthly category
sacut{k}.ye = year(sacut{k}.t);
sacut{k}.mo = month(sacut{k}.t);
sacut{k}.ym=str2num([num2str(sacut{k}.ye),num2str(sacut{k}.mo,'%0.2i')]);
end


% SA data: find group by catagory
% mean by monthly category
for k=1:height(TGinfo)
sa_mean{k}=table();
[x,tm] = findgroups(sacut{k}.ym);

sa_mean{k}.dt=splitapply(@mean, sacut{k}.dt,x);
% sa_mean{k}.dttg=splitapply(@mean, sacut{k}.tgdt,x);
sa_mean{k}.tm=tm;
clearvars x tm

% adding back the time
sa_mean{k}.ye=round((sa_mean{k}.tm)/100,0);
sa_mean{k}.mo=round((((sa_mean{k}.tm)/100)-round((sa_mean{k}.tm)/100,0))*100,0);
sa_mean{k}.time=datetime(sa_mean{k}.ye,sa_mean{k}.mo,1);

%Perform Differences - detrend 
% sa_mean{k}.sadt_Differences=nan(height(sa_mean{k}),1);
% sa_mean{k}.sadt_Differences(2:end)= diff(sa_mean{k}.sadt_Linear);
% me.sadt_Linear=splitapply(@mean, sacut{k}.sadt_Linear,x);
% sa_mean{k}.tg_detrend=detrend(sa_mean{k}.dttg);
sa_mean{k}.sa_detrend=detrend(sa_mean{k}.dt);
end

for k=1:height(TGinfo)
    sa_mean{k}.t=decyear(sa_mean{k}.time);
end


%% 

% TG data: find group by catagory
% for forcasting comparision

%sort by date & monthly category
for k=5:length(tg)
tg{k} = sortrows(tg{k},'Time','ascend');
tg{k}.y = year(tg{k}.Time);
tg{k}.mo = month(tg{k}.Time);
tg{k}.da = day(tg{k}.Time);
tg{k}.ym=str2num([num2str(tg{k}.y),num2str(tg{k}.mo,'%0.2i')]);
tg{k}.ym2=str2num([num2str(tg{k}.y),num2str(tg{k}.mo,'%0.2i'),num2str(tg{k}.da,'%0.2i')]);
end



% monthly mean
for k=1:length(tg)
tg_mean{k}=table();
[x,tm] = findgroups(tg{k}.ym);
tg_mean{k}.dttg=splitapply(@mean, tg{k}.dt2,x); %tide system corrected
tg_mean{k}.tm=tm;
clearvars x tm

% adding back the time
tg_mean{k}.ye=round((tg_mean{k}.tm)/100,0);
tg_mean{k}.mo=round((((tg_mean{k}.tm)/100)-round((tg_mean{k}.tm)/100,0))*100,0);
tg_mean{k}.time=datetime(tg_mean{k}.ye,tg_mean{k}.mo,1);
tg_mean{k}.tg_detrend=detrend(tg_mean{k}.dttg);
end

% daily mean
for k=1:length(tg)
tg_mean_day{k}=table();
[x,tm] = findgroups(tg{k}.ym2);
tg_mean_day{k}.dttg=splitapply(@mean, tg{k}.dt2,x); %tide system corrected
tg_mean_day{k}.tm=tm;
clearvars x tm

tg_mean_day{k}.ye=round((tg_mean_day{k}.tm)/10000,0);
tg_mean_day{k}.mo=round( (((tg_mean_day{k}.tm)/10000)-round((tg_mean_day{k}.tm )/10000,0))*100,0 );
tg_mean_day{k}.da=tg_mean_day{k}.tm-str2num([num2str(tg_mean_day{k}.ye),num2str(tg_mean_day{k}.mo,'%0.2i'),repmat(num2str(00,'%0.2i'),height(tg_mean_day{k}),1)]);
tg_mean_day{k}.time=datetime(tg_mean_day{k}.ye,tg_mean_day{k}.mo,tg_mean_day{k}.da);
tg_mean_day{k}.tg_detrend=detrend(tg_mean_day{k}.dttg);
end

% time to dec year
for k=1:length(tg_mean)
   tg_mean{k}.t=decyear(tg_mean{k}.time); 
   tg_mean_day{k}.t=decyear(tg_mean_day{k}.time); 
end


%% plot TG only
close all
for i=1:6
    figure(i)
plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),tg_mean{i}.dttg(tg_mean{i}.ye>=1995),'-','color',[.7 .7 .7],'DisplayName','TG observation','LineWidth',1.5)
hold on
tr1=fitlm(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),tg_mean{i}.dttg(tg_mean{i}.ye>=1995));
plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),tr1.Fitted,'-k','DisplayName',strcat('TG_{Linear Trend}: ',num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5)

tr2=fitlm(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995&tg_mean{i}.ye<2000)),tg_mean{i}.dttg(tg_mean{i}.ye>=1995&tg_mean{i}.ye<2000));
tr3=fitlm(decyear(tg_mean{i}.time(tg_mean{i}.ye>=2000&tg_mean{i}.ye<2005)),tg_mean{i}.dttg(tg_mean{i}.ye>=2000&tg_mean{i}.ye<2005));
tr4=fitlm(decyear(tg_mean{i}.time(tg_mean{i}.ye>=2005&tg_mean{i}.ye<2010)),tg_mean{i}.dttg(tg_mean{i}.ye>=2005&tg_mean{i}.ye<2010));
tr5=fitlm(decyear(tg_mean{i}.time(tg_mean{i}.ye>=2010&tg_mean{i}.ye<2015)),tg_mean{i}.dttg(tg_mean{i}.ye>=2010&tg_mean{i}.ye<2015));
tr6=fitlm(decyear(tg_mean{i}.time(tg_mean{i}.ye>=2015&tg_mean{i}.ye<2020)),tg_mean{i}.dttg(tg_mean{i}.ye>=2015&tg_mean{i}.ye<2020));


plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995&tg_mean{i}.ye<2000)),tr2.Fitted,'-b','DisplayName',strcat('TG_{95-00}: ',num2str(tr2.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5)
plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=2000&tg_mean{i}.ye<2005)),tr3.Fitted,'-r','DisplayName',strcat('TG_{00-05}: ',num2str(tr3.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5)
plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=2005&tg_mean{i}.ye<2010)),tr4.Fitted,'-g','DisplayName',strcat('TG_{05-10}: ',num2str(tr4.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5)
plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=2010&tg_mean{i}.ye<2015)),tr5.Fitted,'-y','DisplayName',strcat('TG_{10-15}: ',num2str(tr5.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5)
plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=2015&tg_mean{i}.ye<2020)),tr6.Fitted,'-m','DisplayName',strcat('TG_{15-20}: ',num2str(tr6.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5)


xlim([min(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)))-0.1,max(decyear(tg_mean{i}.time(tg_mean{i}.ye<=2022)))+0.1])

legend show

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=20; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',20);

end

close all
for i=1:6
    figure(i)
    
    % all data
    subplot(3,1,1)
    plot(tg{i}.Time,tg{i}.dt2,'.','color',[.7 .7 .7])
    hold on
    tr1=fitlm(decyear(tg{i}.Time),tg{i}.dt2);
    h(1)=plot(tg{i}.Time,tr1.Fitted,'-k','DisplayName',strcat('TG_{Linear Trend}: ',...
        num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5);
    leg = legend(h(1));
    title(leg,'hourly')
    
    %daily mean
    subplot(3,1,2)
    plot(tg_mean_day{i}.time,tg_mean_day{i}.dttg,'.','color',[.7 .7 .7])
    hold on
    tr2=fitlm(decyear(tg_mean_day{i}.time),tg_mean_day{i}.dttg);
    h(1)=plot(tg_mean_day{i}.time,tr2.Fitted,'-k','DisplayName',strcat('TG_{Linear Trend}: ',...
        num2str(tr2.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5);
    leg = legend(h(1));
    title(leg,'daily')
    
    %monthly mean
    subplot(3,1,3)
    plot(tg_mean{i}.time,tg_mean{i}.dttg,'.','color',[.7 .7 .7])
    hold on
    tr3=fitlm(decyear(tg_mean{i}.time),tg_mean{i}.dttg);
    h(1)=plot(tg_mean{i}.time,tr3.Fitted,'-k','DisplayName',strcat('TG_{Linear Trend}: ',....
        num2str(tr3.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5);
    leg = legend(h(1));
    title(leg,'monthly')   
        
        
end

%% load data
load('F:\Analysis\BS\trend\ASL trend\data_4st_allsa.mat')

%% 2nd deg poly

close all
for i=1:4
figure(i)
h(1)=plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),tg_mean{i}.dttg(tg_mean{i}.ye>=1995),'-','color',[.7 .7 .7],'DisplayName','TG observation','LineWidth',1.5);
hold on
h(2)=plot(decyear(sa_mean{i}.time(sa_mean{i}.ye>=1995)),sa_mean{i}.dt(sa_mean{i}.ye>=1995),'-','color',[0.3010 0.7450 0.9330],'DisplayName','SA observation','LineWidth',1.5);

tr1=fitlm(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),tg_mean{i}.dttg(tg_mean{i}.ye>=1995));


[p1,s1] = polyfit(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),tg_mean{i}.dttg(tg_mean{i}.ye>=1995),1);
[y_fit1,f1]= polyval(p1,decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),s1);
h(3)=plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),y_fit1,'-','color',[0.4940 0.1840 0.5560],'LineWidth',2.5,'DisplayName',...
    strcat('Linear Trend= ',num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'));
%     strcat('TG_{Linear Trend}: ',num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'));

plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),y_fit1+f1,'--','color',[0.4940 0.1840 0.5560])
plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),y_fit1-f1,'--','color',[0.4940 0.1840 0.5560])


[p2,s2] = polyfit(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),tg_mean{i}.dttg(tg_mean{i}.ye>=1995),2);
[y_fit2,f2]= polyval(p2,decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),s2);
h(4)=plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),y_fit2,'-g','LineWidth',2.5,'DisplayName',strcat('2^{nd}poly trend'));
% plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),y_fit2+f2,'k--',decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),y_fit2-f2,'k--')
% 


[p3,s3] = polyfit(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),tg_mean{i}.dttg(tg_mean{i}.ye>=1995),5);
[y_fit3,f3]= polyval(p3,decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),s3);
h(5)=plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),y_fit3,'-k','LineWidth',2.5,'DisplayName',strcat('5^{th}poly trend'));
% plot(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),y_fit2+f2,'k--',decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)),y_fit2-f2,'k--')
% 


xlim([min(decyear(tg_mean{i}.time(tg_mean{i}.ye>=1995)))-0.1,max(decyear(tg_mean{i}.time(tg_mean{i}.ye<=2022)))+0.1])

legend(h(1:5))

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=20; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',20);
end

%% 

% close all

for i=11:14
% for i=1:4
figure(i)

mintime=min(sa_mean{i}.t);
maxtime=max(sa_mean{i}.t);
ti=sa_mean{i}.t;
dttg=interp1(tg_mean{i}.t,tg_mean{i}.dttg,ti);
dtsa=sa_mean{i}.dt;


tr1=fitlm(ti,dttg);
tr2=fitlm(ti,dtsa);


h(1)=plot(ti,tr1.Fitted,'-k','LineWidth',2.5,'DisplayName',...
    strcat('TG_{Linear Trend}= ',num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'));

hold on
h(2)=plot(ti,tr2.Fitted,'-','color',[0.4940 0.1840 0.5560],'LineWidth',2.5,'DisplayName',...
    strcat('SA_{Linear Trend}= ',num2str(tr2.Coefficients.Estimate(2)*10,2),' mm/year'));

% plot(ti,y_fit1+delta,'--','color',[0.4940 0.1840 0.5560])
% plot(ti,y_fit1-delta,'--','color',[0.4940 0.1840 0.5560])


p1= polyfit(ti,dttg,2);
fit1=polyval(p1,ti);
h(3)=plot(ti,fit1,'-','color',[0 0.4470 0.7410],'LineWidth',2.5,'DisplayName',strcat('2^{nd}poly trend'));

p2= polyfit(ti,dttg,3);
fit2=polyval(p2,ti);
h(4)=plot(ti,fit2,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2.5,'DisplayName',strcat('3^{nd}poly trend'));

p3= polyfit(ti,dttg,4);
fit3=polyval(p3,ti);
h(5)=plot(ti,fit3,'-','color',[0.9290 0.6940 0.1250],'LineWidth',2.5,'DisplayName',strcat('4^{th}poly trend'));

p4= polyfit(ti,dttg,5);
fit4=polyval(p4,ti);
h(6)=plot(ti,fit4,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2.5,'DisplayName',strcat('5^{th}poly trend'));

p5= polyfit(ti,dttg,6);
fit5=polyval(p5,ti);
h(7)=plot(ti,fit5,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2.5,'DisplayName',strcat('6^{th}poly trend'));



xlim([mintime-0.5,maxtime+0.5])

legend(h(1:7),'Location','southeast')

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=20; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
xticks(1995:1:2020)
set(gca,'fontname','Times New Roman','FontSize',20);

clearvars mintime maxtime ti dtsa dttg h p1 p2 p3 p4 p5 fit1 fit2 fit3 fit4 fit5 tr1 tr2 delta
end