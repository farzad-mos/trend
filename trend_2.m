cd F:\Analysis\BS\trend

clear
clc

load('TG.mat')
% load('F:\Analysis\BS\trend\data\sa\all\SA_orig.mat')
load('sa.mat')


%1   2   3   4   5   6   7   8   9           
%CS2 ENV ER2 JA1 JA2 JA3 SRL S3A S3B


%%  Geoplot study area

geoplot(TGinfo.Lat,TGinfo.Lon,'^k','MarkerFaceColor','k','MarkerSize',10)
hold on
% geoplot(sa.lat(sa.said==2&sa.cycle==20),sa.lon(sa.said==2&sa.cycle==20),'.','color',[0.4660 0.6740 0.1880])
% geoplot(sa.lat(sa.said==4&sa.cycle==17),sa.lon(sa.said==4&sa.cycle==17),'.','color',[0 0.4470 0.7410]) 

geoplot(sa.lat(sa.said==7),sa.lon(sa.said==7),'.','color',[0.4660 0.6740 0.1880])
geoplot(sa.lat(sa.said==8),sa.lon(sa.said==8),'.','color',[0 0.4470 0.7410]) 

for i=1:height(TGinfo)
    [lat,lon] = scircle1(TGinfo.Lat(i),TGinfo.Lon(i),km2deg(50),[],[],[],50);
geoplot(lat,lon,'r')
clearvars lat lon
end

ax=gca; ax.FontSize=14; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');
clearvars ax i 

%% remove gross errors & outliers 

%convert SSH to DT
sa.dt=sa.ssh-sa.N+sa.dac;


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
%% SA corrections
% tide system conversion (zero tide to mean tide) Ekman 1989

% ******** apply once *********
sa.sadt=sa.dt+(0.09940-0.29541*((sind(sa.lat)).^2)-0.00042*sind(sa.lat).^4); %******* check ****


% Ellipsoidal height correction + TRF conversion
% load T/P ellipsoid
load('TP.mat')
spheroid = referenceEllipsoid('GRS 80');

% sa.sadt=nan(height(sa),1);

for ye=min(sa.year_sa):max(sa.year_sa)

    % convert Lat, Lon, H to X Y Z (from T/P Ellipsoid)
    [x,y,z] = geodetic2ecef(TP,sa.lat(sa.year_sa==ye),sa.lon(sa.year_sa==ye),sa.sadt(sa.year_sa==ye));
    
    % transfering form: ITRF2008	to: ETRF2000	 epoch: SA Year;
    X2=itrstrafo([x,y,z],'ITRF2008',ye,'ETRF2000',ye);
    
    % convert X Y Z to Lat Lon H (to GRS 80 Ellipsoid)
    [~,~,h] = ecef2geodetic(spheroid,X2(:,1),X2(:,2),X2(:,3));
        
    sa.sadt(sa.year_sa==ye)=h*100;
    clearvars h X2 x y z
end



% elipsoid correction
[x,y,z] = geodetic2ecef(TP,sa.lat,sa.lon,sa.sadt);
[~,~,h] = ecef2geodetic(spheroid,x,y,z);
sa.sadt2(:)=h*100;

sa(abs(sa.sadt)>=200,:)=[];
    
clearvars ye X2 x y z


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

%% combine all TG tables

for i=1:height(TGinfo)
T=eval(['TG' num2str(i)]);
[~,ia] = unique(T.date);
tg{i}=T(ia,:);
clearvars T ia
end
clear i
%% mean SA at eact TG per cycle
for k=1:height(TGinfo)
    sacor=table();
    H=1;
    for i=min(sa.said):max(sa.said)
        
        id=unique(sa_tg{1, k}.id(sa_tg{1, k}.said==i));
        if ~isempty(id)
            for j=id(1):id(end)
                
                if ~isempty(sa_tg{1, k}.ssh(sa_tg{1, k}.said==i&sa_tg{1, k}.id==j))
                    
                    sacor.dt(H)=mean((sa_tg{1, k}.sadt(sa_tg{1, k}.said==i&sa_tg{1, k}.id==j)),'omitnan');
                    sacor.dac(H)=mean((sa_tg{1, k}.dac(sa_tg{1, k}.said==i&sa_tg{1, k}.id==j)),'omitnan');
                    time=sa_tg{1, k}.t(sa_tg{1, k}.said==i&sa_tg{1, k}.id==j);
                    sacor.t(H)=time(1);
                    sacor.said(H)=i;
                    sacor.tgdt(H)=interp1(tg{1, k}.Time,tg{1, k}.dt,time(1));
                    H=H+1;
                    clear time
                end
            end
        end
    end
    satg{k}=sacor;
    clear sacor
end

%% bias removal
% remove +- 100 cm dt
for k=1:height(TGinfo)
satg{1,k}((abs(satg{1, k}.dt)>=100),:)=[];
end

% apply mean diff
for k=1:height(TGinfo)
satg{1, k}.sadt=(satg{1, k}.dt)+(mean(satg{1, k}.tgdt,'omitnan')-mean(satg{1, k}.dt,'omitnan'));
end

for k=1:height(TGinfo)
rmse(k)=rms(((satg{1, k}.tgdt)-satg{1, k}.sadt),'omitnan');
end
%% plot sa tg
% first 5 TGs
% only 5 SA ***** change *****

for k=1:5
    
    subplot(5,1,k)
    plot(tg{1, k}.Time,tg{1, k}.dt,'-','color',[.7 .7 .7])
    hold on
    plot(satg{1, k}.t(satg{1, k}.said==2),satg{1, k}.sadt(satg{1, k}.said==2),'*','color',[0.4940 0.1840 0.5560])
    plot(satg{1, k}.t(satg{1, k}.said==3),satg{1, k}.sadt(satg{1, k}.said==3),'*','color',[0.4660 0.6740 0.1880])
    plot(satg{1, k}.t(satg{1, k}.said==7),satg{1, k}.sadt(satg{1, k}.said==7),'*','color',[0.6350 0.0780 0.1840])
    plot(satg{1, k}.t(satg{1, k}.said==8),satg{1, k}.sadt(satg{1, k}.said==8),'*','color',[0 0.4470 0.7410])
    plot(satg{1, k}.t(satg{1, k}.said==9),satg{1, k}.sadt(satg{1, k}.said==9),'*','color',[0.1 0.8 1])

    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
    ylabel('DT [cm]','FontSize',18,'FontWeight','bold');
    ylim([-100,200])
    if k~=5
        xticks([])
    end
    xlim([datetime(1995,1,1),datetime(2020,1,1)])
    
end

clearvars ax k 


%% trend

for k=1:5
    p1 = polyfit(tg{1, k}.date,tg{1, k}.dt,1);
    f1 = polyval(p1,tg{1, k}.date);
    
    p2 = polyfit(decyear(satg{1, k}.t),satg{1, k}.sadt,1);
    f2 = polyval(p2,decyear(satg{1, k}.t));
    
    p3 = polyfit(decyear(satg{1, k}.t),satg{1, k}.sadt-satg{1, k}.dac,1);
    f3 = polyval(p3,decyear(satg{1, k}.t));
    
    tr1=fitlm(tg{1, k}.date,tg{1, k}.dt);
    tr2=fitlm(decyear(satg{1, k}.t),satg{1, k}.sadt);
    tr3=fitlm(decyear(satg{1, k}.t),satg{1, k}.sadt-satg{1, k}.dac);
    
    
    subplot(5,1,k)
    % plot(tg{1, k}.date,tg{1, k}.dt,'.','color',[.9 .9 .9]);
    % hold on
    % plot(satg{1, k}.t,satg{1, k}.sadt,'.r');
    
    p(1)=plot(tg{1, k}.Time,f1,'-k','LineWidth',2.5,'DisplayName',strcat('TG trend= ',num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'));
    hold on
    p(2)=plot(satg{1, k}.t,f2,'-r','LineWidth',2.5,'DisplayName',strcat('SA trend= ',num2str(tr2.Coefficients.Estimate(2)*10,2),' mm/year'));
    p(3)=plot(satg{1, k}.t,f3,'-b','LineWidth',2.5,'DisplayName',strcat('SA_D_A_C trend= ',num2str(tr3.Coefficients.Estimate(2)*10,2),' mm/year'));
    
    
    
    legend(p(1:3),'FontSize',10,'Location','northwest')
    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
    ylabel('DT [cm]','FontSize',18,'FontWeight','bold');
    
    clearvars p ax f1 f2 f3 p1 p2 p3 tr1 tr2 tr3
    
    if k~=5
        xticks([])
    end
end
clear k

%% select TG data
clear res
clc

% select tg and associated data
k=1;


res=table();
H=1;
dt_sa=satg{1, k}.sadt;
dt_tg=tg{1, k}.dt;
date_tg=tg{1, k}.date;



%% simulate data

% Autoregression (AR)

for ar=10:5:30 % 10, 15, 20, 25, 30
sys = arima(ar,0,0);
Md1 = estimate(sys,dt_sa);
% summarize(Md1)

residual1 = infer(Md1,dt_tg);
prediction1 = dt_tg + residual1;

R=corrcoef(dt_tg,prediction1);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction1);
res.RMSE(H)=rms(dt_tg-prediction1);
res.model(H,:)=string(strcat('AR(',num2str(ar),')'));
H=H+1;
end

% Moving Average
for  ma=10:5:25  %10 15 20 25
sys = arima(0,0,ma);
Md1 = estimate(sys,dt_sa);
residual2 = infer(Md1,dt_tg);
prediction2 = dt_tg + residual2;

R = corrcoef(dt_tg,prediction2);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction2);
res.RMSE(H)=rms(dt_tg-prediction2);
res.model(H,:)=string(strcat('MA(',num2str(ma),')'));

H=H+1;
end

% Autoregressive Moving Average (ARMA)

for ar=10:5:25 %10, 15, 20, 25
ma=2; %2, 5

sys = arima(ar,0,ma);
Md1 = estimate(sys,dt_sa);
residual3 = infer(Md1,dt_tg);
prediction3 = dt_tg + residual3;

R = corrcoef(dt_tg,prediction3);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction3);
res.RMSE(H)=rms(dt_tg-prediction3);
res.model(H,:)=string(strcat('ARMA(',num2str(ar),',',num2str(ma),')'));
H=H+1;
end

for ar=10:5:25 %10, 15, 20, 25
ma=5; %2, 5

sys = arima(ar,0,ma);
Md1 = estimate(sys,dt_sa);
residual3 = infer(Md1,dt_tg);
prediction3 = dt_tg + residual3;

R = corrcoef(dt_tg,prediction3);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction3);
res.RMSE(H)=rms(dt_tg-prediction3);
res.model(H,:)=string(strcat('ARMA(',num2str(ar),',',num2str(ma),')'));
H=H+1;
end

% Autoregressive Integrated Moving Average (ARIMA)
for ar=10:5:25 %10, 15, 20, 25

ma=2; %2, 5
I=5; 
sys = arima(ar,I,ma);
Md1 = estimate(sys,dt_sa);
residual4 = infer(Md1,dt_tg);
prediction4 = dt_tg + residual4;

R = corrcoef(dt_tg,prediction4);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction4);
res.RMSE(H)=rms(dt_tg-prediction4);
res.model(H,:)=string(strcat('ARIMA(',num2str(ar),',',num2str(I),',',num2str(ma),')'));
H=H+1;
end

for ar=10:5:25 %10, 15, 20, 25

ma=5; %2, 5
I=5; 
sys = arima(ar,I,ma);
Md1 = estimate(sys,dt_sa);
residual4 = infer(Md1,dt_tg);
prediction4 = dt_tg + residual4;

R = corrcoef(dt_tg,prediction4);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction4);
res.RMSE(H)=rms(dt_tg-prediction4);
res.model(H,:)=string(strcat('ARIMA(',num2str(ar),',',num2str(I),',',num2str(ma),')'));
H=H+1;
end

for S=12:12:24
% Seasonal Autoregressive Integrated Moving-Average (SARIMA)
sys = arima('Constant',NaN,'ARLags',1:4,'D',0,'MALags',1:2,'SARLags',[12,24,36,48],'Seasonality',S,'SMALags',S,'Distribution','Gaussian');
Md1 = estimate(sys,dt_sa);
residual5 = infer(Md1,dt_tg);
prediction5 = dt_tg + residual5;


R = corrcoef(dt_tg,prediction5);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction5);
res.RMSE(H)=rms(dt_tg-prediction5);
res.model(H,:)=string(strcat('SARIMA(',num2str(S),')'));
H=H+1;
end

% GARCH Model
GARCH_X1 = garch('Offset',0,'GARCHLags',1:14,'ARCHLags',50,'Distribution','Gaussian');
GARCH_X1 = estimate(GARCH_X1,dt_sa,'Display','off');
residual6 = infer(GARCH_X1,dt_tg);
prediction6 = dt_tg + residual6;

R = corrcoef(dt_tg,prediction6);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction6);
res.RMSE(H)=rms(dt_tg-prediction6);
res.model(H,:)="GARCH";
H=H+1;


%Glostan, Jagannathan and Runkle GARCH Model
% this give the same result as GRACH model in resampling
GJR_X1 = gjr('Offset',0,'GARCHLags',1:3,'ARCHLags',1,'LeverageLags',1,'Distribution','Gaussian');
GJR_X1 = estimate(GJR_X1,dt_sa,'Display','off');
residual7 = infer(GJR_X1,dt_tg);
prediction7 = dt_tg + residual7;

R = corrcoef(dt_tg,prediction7);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction7);
res.RMSE(H)=rms(dt_tg-prediction7);
res.model(H,:)="GJR";
H=H+1;



%% select best case
for k=1:height(TGinfo)

    dt_sa=satg{1, k}.sadt;
    dt_tg=tg{1, k}.dt;
    date_tg=tg{1, k}.date;
    
% Autoregression (AR)    
%     ar= 20;
%     sys = arima(ar,0,0);
%     Md1 = estimate(sys,dt_sa);
%     residual1 = infer(Md1,dt_tg);
%     prediction1 = dt_tg + residual1;

% SARIMA
S=24;
sys = arima('Constant',NaN,'ARLags',1:4,'D',0,'MALags',1:2,'SARLags',[12,24,36,48],'Seasonality',S,'SMALags',S,'Distribution','Gaussian');
Md1 = estimate(sys,dt_sa);
residual5 = infer(Md1,dt_tg);
prediction1 = dt_tg + residual5;
        
    subplot(4,4,k)
    plot(date_tg,dt_tg-prediction1,'color',[0.3010 0.7450 0.9330])
    yline(0,'-','Refrence','fontname','Times New Roman','FontSize',13,'FontWeight','Bold','LineWidth',2);
   if k==1 || k==5 || k==9
    ylabel('\DeltaDT [cm]','FontSize',18,'FontWeight','bold');
    end
    ylim([-100 50])
    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
    set(gca,'fontname','Times New Roman','FontSize',18);
    
    xlim([1995 2020])
    tg{1, k}.dtpredict=prediction1;
    clear prediction1   
   

end

%% box plot
    
for k=1:height(TGinfo)
    meanWeight = groupsummary(tg{1, k}.dt-tg{1, k}.dtpredict,month(tg{1, k}.Time),'mean');
    
    subplot(3,4,k)
    boxchart(month(tg{1, k}.Time),tg{1, k}.dt-tg{1, k}.dtpredict,'MarkerStyle','none','WhiskerLineColor','none')
    hold on
    
    plot(meanWeight,'-b')
    yline(0,'--k','LineWidth',1.5);
    ylim([-15 15])
    xlim([0 13])
    
    if k==1 || k== 5 || k==9
        xlabel('Month','FontSize',18,'FontWeight','bold');
        ylabel('\DeltaDT (TG-SA_r_e_s_a_m_p_l_e_d) [cm]','FontSize',18,'FontWeight','bold');
    end
    xticks(1:1:12)
    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
    set(gca,'fontname','Times New Roman','FontSize',18);
    box on
    
    mae(tg{1, k}.dt-tg{1, k}.dtpredict)
    rms(tg{1, k}.dt-tg{1, k}.dtpredict)
    
    if k~=3
        xticks([])
    end
    
end
    
%% trend

for k=1:height(TGinfo)
    p1 = polyfit(tg{1, k}.date,tg{1, k}.dt,1);
    f1 = polyval(p1,tg{1, k}.date);
    
    p2 = polyfit(decyear(satg{1, k}.t),satg{1, k}.sadt-satg{1, k}.dac,1);
    f2 = polyval(p2,decyear(satg{1, k}.t));
    
    p3 = polyfit(tg{1, k}.date,tg{1, k}.dtpredict,1);
    f3 = polyval(p3,tg{1, k}.date);
    
    tr1=fitlm(tg{1, k}.date,tg{1, k}.dt);
    tr2=fitlm(decyear(satg{1, k}.t),satg{1, k}.sadt-satg{1, k}.dac);
    tr3=fitlm(tg{1, k}.date,tg{1, k}.dtpredict);
    
    subplot(3,4,k)
    % plot(tg{1, k}.date,tg{1, k}.dt,'.','color',[.9 .9 .9]);
    % hold on
    % plot(satg{1, k}.t,satg{1, k}.sadt,'.r');
    
    p(1)=plot(tg{1, k}.date,f1,'-k','LineWidth',2.5,'DisplayName',strcat('TG trend= ',num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'));
    hold on
    p(2)=plot(decyear(satg{1, k}.t),f2,'-','color',[0.3010 0.7450 0.9330],'LineWidth',2.5,'DisplayName',strcat('SA trend= ',num2str(tr2.Coefficients.Estimate(2)*10,2),' mm/year'));
    p(3)=plot(tg{1, k}.date,f3,'-b','LineWidth',2.5,'DisplayName',strcat('SA_r_e_s_a_m_p_l_e_d trend= ',num2str(tr3.Coefficients.Estimate(2)*10,2),' mm/year'));
    
    
    legend(p(1:3),'FontSize',12,'Location','southeast')
    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  %grid on;
    
    if k==1 || k==5 || k==9
        ylabel('DT [cm]','FontSize',18,'FontWeight','bold');
    end
    
    if k==1 || k==2 || k==3 || k==4 || k==5 || k==6 || k==7 || k==8 
        xticks([])
    end
        xlim([1995 2020])
        ylim([10 30])

end
%% 
latlim = [53 67];
lonlim = [10 30.5];

ax = usamap(latlim, lonlim);

hold on

geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
setm(gca,'FLineWidth',5,'Grid','on','FontSize',18,'fontweight','bold')
scatterm(TGinfo.Lat,TGinfo.Lon,800,TGinfo.tr,'filled','o','MarkerEdgeColor','k');
c=colorbar; 
c.Label.String = 'Sea Level Trend [mm/year]';
colormap(flipud(hot))
caxis([3 5]); %fix the bar: SA

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; 


% xLoc =1.3001e+05;
% yLoc =7.0876e+06;
% [xLoc, yLoc] = ginput(1);


scaleruler('Units', 'km', 'RulerStyle', 'patches', ...
'XLoc', xLoc, 'YLoc', yLoc,"FontSize",10);
    
%% 
hold on
subplot (2,3,1)
histogram(dt_tg-prediction1,'Normalization','probability')
xlim([-15 15])
subplot (2,3,2)
histogram(dt_tg-prediction2,'Normalization','probability')
xlim([-15 15])
subplot (2,3,3)
histogram(dt_tg-prediction3,'Normalization','probability')
xlim([-15 15])
subplot (2,3,4)
histogram(dt_tg-prediction4,'Normalization','probability')
xlim([-15 15])
subplot (2,3,5)
histogram(dt_tg-prediction5,'Normalization','probability')
xlim([-15 15])
subplot (2,3,6)
histogram(dt_tg-prediction6,'Normalization','probability')
xlim([-15 15])

%% Long Short-Term Memory Networks (LSTM)

numTimeStepsTrain = floor(0.9*numel(dt_sa));
dataTrain = dt_sa(1:numTimeStepsTrain+1);
dataTest = dt_sa(numTimeStepsTrain+1:end);

% Standardize Data

mu = mean(dataTrain);
sig = std(dataTrain);
dataTrainStandardized = (dataTrain - mu) / sig;
% [dataTrainStandardized,inputps] = mapminmax(dataTrain); %2nd option

% Prepare Predictors and Responses
XTrain = dataTrainStandardized(1:end-1)';
YTrain = dataTrainStandardized(2:end)';

% Define LSTM Network Architecture
numFeatures = 1;
numResponses = 1;
numHiddenUnits = 200;

layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits)
    fullyConnectedLayer(numResponses)
    regressionLayer];

% layers2 = [ ...
%     sequenceInputLayer(numFeatures)
%     bilstmLayer(numHiddenUnits)
%     fullyConnectedLayer(numResponses)
%     regressionLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',250, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.005, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.2, ...
    'Verbose',0, ...
    'Plots','training-progress');
% Train LSTM Network
net = trainNetwork(XTrain,YTrain,layers,options);

% Forecast Future Time Steps
dataTestStandardized = (dataTest - mu) / sig;
% dataTestStandardized = mapminmax('apply',dataTest,inputps); %2nd option



XTest = dataTestStandardized(1:end-1);
net = predictAndUpdateState(net,XTrain);
[net,YPred] = predictAndUpdateState(net,YTrain(end));

numTimeStepsTest = numel(XTest);
for i = 2:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,YPred(:,i-1),'ExecutionEnvironment','cpu');
end

% Unstandardize the predictions using the parameters calculated earlier.
YPred = sig*YPred + mu;
% YPred = mapminmax('reverse',YPred,outputps); %2nd option

% Goodness of the fit statistics 
YTest = dataTest(2:end);

RMSE=rms(YPred-YTest')
MAE=mae(YPred-YTest')
R = corrcoef(YPred,YTest')



figure
plot(dataTrain(1:end-1))
hold on
idx = numTimeStepsTrain:(numTimeStepsTrain+numTimeStepsTest);
plot(idx,[dt_sa(numTimeStepsTrain) YPred],'.-')
hold off
ylabel("SLA")
title("Forecast")
legend(["Observed" "Forecast"])

figure
subplot(2,1,1)
plot(YTest)
hold on
plot(YPred,'.-')
hold off
legend(["Observed" "Forecast"])
ylabel("Cases")
title("Forecast")

subplot(2,1,2)
stem(YPred - YTest')
xlabel("Month")
ylabel("Error")
title("RMSE = " + rmse)



% Update Network State with Observed Values
net = resetState(net);
net = predictAndUpdateState(net,XTrain);

YPred = [];
XTest=XTest';
numTimeStepsTest = numel(XTest);
for i = 1:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,XTest(:,i),'ExecutionEnvironment','cpu');
end

% Unstandardize the predictions using the parameters calculated earlier.
YPred = sig*YPred + mu;
% YPred = mapminmax('reverse',YPred,outputps); %2nd option


RMSE=rms(YPred-YTest')
MAE=mae(YPred-YTest')
R = corrcoef(YPred,YTest')

figure
subplot(2,1,1)
plot(YTest)
hold on
plot(YPred,'.-')
hold off
legend(["Observed" "Predicted"])
ylabel("Cases")
title("Forecast with Updates")

subplot(2,1,2)
stem(YPred - YTest')
xlabel("Month")
ylabel("Error")
title("RMSE = " + rmse)

%% Extreme Learning Machine (ELM)
 
 % II. Training set/test set generation
 numTimeStepsTrain= floor(0.9*numel(dt_sa));
 % training set
P_train = dt_sa(1:numTimeStepsTrain+1)';
T_train = dt_sa(1:numTimeStepsTrain+1)';
 
 % test set
P_test = dt_sa(numTimeStepsTrain+1:end)';
T_test = dt_sa(numTimeStepsTrain+1:end)';
N = size(P_test,2);
 
% III. Data normalization
% 1. Training set
[Pn_train,inputps] = mapminmax(P_train);
Pn_test = mapminmax('apply',P_test,inputps);
% 2. Test set
[Tn_train,outputps] = mapminmax(T_train);
Tn_test = mapminmax('apply',T_test,outputps);
 
% IV. ELM creation/training
[IW,B,LW,TF,TYPE] = elmtrain(Pn_train,Tn_train,30,'sig',0);
 
% V. ELM simulation test
tn_sim = elmpredict(Pn_test,IW,B,LW,TF,TYPE);
% 1. Anti-normalization
T_sim = mapminmax('reverse',tn_sim,outputps);
 
% VI. Comparison of results
result = [T_test' T_sim'];
N = length(T_test);

RMSE=rms(T_sim-T_test)
MAE=mae(T_sim-T_test)
MSE=mse(T_sim - T_test);
R = corrcoef(T_sim,T_test)


% VII. Drawing
figure(1)
plot(1:N,T_test,'r-*',1:N,T_sim,'b:o')
grid on
legend('true value', 'predicted value')
Xlabel('sample number')
Ylabel('octane number')
String = {'ELM';['(mse = ' num2str(MSE) ' R^2 = ' num2str(R(1,2)) ')']};
title(string)



