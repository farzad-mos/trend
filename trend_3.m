load('satg.mat')
load('tg1.mat')

% cut data 
for k=1:16
sacut{k}=satg{k}(satg{k}.said<=3,:);
end

%sort by date
for k=1:16
sacut{k} = sortrows(sacut{k},'t','ascend');
end

%monthly category
for k=1:16
sacut{k}.ye = year(sacut{k}.t);
sacut{k}.mo = month(sacut{k}.t);
sacut{k}.ym=str2num([num2str(sacut{k}.ye),num2str(sacut{k}.mo,'%0.2i')]);
end
%% Transform Data into Stationary Data

for k=1:size(sacut,2)
% Log
sacut{k}.sadt_Log = log(sacut{k}.sadt);
% Remove Linear trend
% save the trend to retrend back the data
sacut{k}.sadt_logLinear=detrend(sacut{k}.sadt_Log);
sacut{k}.sadt_Linear=detrend(sacut{k}.sadt);
%Perform Differences
sacut{k}.sadt_Differences=nan(height(sacut{k}),1);
sacut{k}.sadt_Differences(2:end)= diff(sacut{k}.sadt_Linear);
end

%% Test Stationary
close all

for k=1:14

% Hypothesis test : Augmented Dickey-filer, KPSS, Leybourne-McCabe, Philip-Peron, Variance Ratio
% Hypothesis test : Egale's ARCH, Ljung-Box Q-test
[h,p] = adftest(sacut{k}.sadt_Differences);
% if test is valid, it will return 1 or else 0
% P-Value will also verify the stationarity of data, It may get extremely low if your data is not valid.
display(h)


% Autocorrelaton function (ACF)
% Identify series with serial correlation
% Determine whether an AR model is apprpriate
% Identify significant MA lags for model identification
figure(k)
subplot(3,3,1)
plot(sacut{k}.t,sacut{k}.sadt,'-','color',[0.4940 0.1840 0.5560])
hold on
p1 = polyfit(decyear(sacut{k}.t),sacut{k}.sadt,1);
f1 = polyval(p1,decyear(sacut{k}.t));
tr1=fitlm(decyear(sacut{k}.t),sacut{k}.sadt);
p(1)=plot(sacut{k}.t,f1,'-k','LineWidth',2.5,'DisplayName',strcat('trend= ',num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'));
h=yline(mean(sacut{k}.sadt),'--k',strcat('MDT=',num2str(mean(sacut{k}.sadt),2),' [cm]'),'LineWidth',2);
h.LabelHorizontalAlignment='left';
h.FontSize=18;
legend(p(1),'FontSize',18,'Location','northwest')
legend boxoff
ylabel('DT [cm]')
xlabel([])
title('DT')

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);

subplot(3,3,2)
autocorr(sacut{k}.sadt,'NumLags',20,'NumMA',1)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);
xlabel([])
ylabel('ACF')
title('DT')

% Partial ACF (PACF)
% Identify Series with Serial Correlation
% Determine whether an MA model is appropriate
% Identify significant AR lags for model identification.
subplot(3,3,3)
parcorr(sacut{k}.sadt,'NumAR',1,'NumLags',20)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);
xlabel([])
ylabel('PACF')
title('DT')

subplot(3,3,4)
plot(sacut{k}.t,sacut{k}.sadt_Linear,'-','color',[0.4940 0.1840 0.5560])
hold on
h=yline(mean(sacut{k}.sadt_Linear),'--k',strcat('Mean detrended=',num2str(mean(sacut{k}.sadt_Linear),2),' [cm]'),'LineWidth',2);
h.LabelHorizontalAlignment='left';
h.FontSize=18;
ylabel('DT [cm]')
xlabel([])
title('detrended DT')
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);

%linear
subplot(3,3,5)
autocorr(sacut{k}.sadt_Linear,'NumLags',20,'NumMA',1)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);
xlabel([])
ylabel('ACF')
title('detrended DT')

subplot(3,3,6)
parcorr(sacut{k}.sadt_Linear,'NumAR',1,'NumLags',20)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);
xlabel([])
ylabel('PACF')
title('detrended DT')


subplot(3,3,7)
plot(sacut{k}.t,sacut{k}.sadt_Differences,'-','color',[0.4940 0.1840 0.5560])
hold on
h=yline(mean(sacut{k}.sadt_Differences,'omitnan'),'--k',strcat('Mean diff=',num2str(mean(sacut{k}.sadt_Differences,'omitnan'),2),' [cm]'),'LineWidth',2);
h.LabelHorizontalAlignment='left';
h.FontSize=18;
ylabel('DT [cm]')
xlabel('date')
title('diff detrended DT')
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);


%difference
subplot(3,3,8)
autocorr(sacut{k}.sadt_Differences,'NumLags',20,'NumMA',1)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);
ylabel('ACF')
title('diff detrended DT')

subplot(3,3,9)
parcorr(sacut{k}.sadt_Differences,'NumAR',1,'NumLags',20)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);
ylabel('PACF')
title('diff detrended DT')

end

%% save the results

% select TG data
clear res
clc

% select tg and associated data
k=4;


res=table();
H=1;
dt_sa=sacut{k}.sadt_Linear;
dt_tg=tg{1, k}.dt;
date_tg=tg{1, k}.date;



%% simulate data

% Autoregression (AR)

for ar=1:1:10 % 10, 15, 20, 25, 30
sys = arima(ar,0,0);
[Md1 ,~,Loglikehood] = estimate(sys,dt_sa);


% summarize(Md1)

residual1 = infer(Md1,dt_tg);
prediction1 = dt_tg + residual1;

R=corrcoef(dt_tg,prediction1);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction1);
res.RMSE(H)=rms(dt_tg-prediction1);
res.model(H,:)=string(strcat('AR(',num2str(ar),')'));

% Goodness-of-Fit Checks
% Alaike or Bayesoan: to avoid the fitting is over-fit.
[a,b] = aicbic(Loglikehood,2,250);
res.aic(H,:)=a;
res.bic(H,:)=b;
clearvars a b Loglikehood Md1 sys R residual1 prediction1 


H=H+1;
end

% Moving Average
for  ma=5:5:25  %10 15 20 25
sys = arima(0,0,ma);
[Md1 ,~,Loglikehood] = estimate(sys,dt_sa);
residual2 = infer(Md1,dt_tg);
prediction2 = dt_tg + residual2;

R = corrcoef(dt_tg,prediction2);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction2);
res.RMSE(H)=rms(dt_tg-prediction2);
res.model(H,:)=string(strcat('MA(',num2str(ma),')'));


% Goodness-of-Fit Checks
% Alaike or Bayesoan: to avoid the fitting is over-fit.
[a,b] = aicbic(Loglikehood,2,250);
res.aic(H,:)=a;
res.bic(H,:)=b;
clearvars a b Loglikehood Md1 sys R residual2 prediction2 
H=H+1;
end

% Autoregressive Moving Average (ARMA)

for ar=5:5:25 %10, 15, 20, 25
ma=2; %2, 5

sys = arima(ar,0,ma);
[Md1 ,~,Loglikehood] = estimate(sys,dt_sa);
residual3 = infer(Md1,dt_tg);
prediction3 = dt_tg + residual3;

R = corrcoef(dt_tg,prediction3);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction3);
res.RMSE(H)=rms(dt_tg-prediction3);
res.model(H,:)=string(strcat('ARMA(',num2str(ar),',',num2str(ma),')'));

% Goodness-of-Fit Checks
% Alaike or Bayesoan: to avoid the fitting is over-fit.
[a,b] = aicbic(Loglikehood,2,250);
res.aic(H,:)=a;
res.bic(H,:)=b;
clearvars a b Loglikehood Md1 sys R residual3 prediction3 

H=H+1;
end

for ar=5:5:25 %10, 15, 20, 25
ma=5; %2, 5

sys = arima(ar,0,ma);
[Md1 ,~,Loglikehood] = estimate(sys,dt_sa);
residual3 = infer(Md1,dt_tg);
prediction3 = dt_tg + residual3;

R = corrcoef(dt_tg,prediction3);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction3);
res.RMSE(H)=rms(dt_tg-prediction3);
res.model(H,:)=string(strcat('ARMA(',num2str(ar),',',num2str(ma),')'));

% Goodness-of-Fit Checks
% Alaike or Bayesoan: to avoid the fitting is over-fit.
[a,b] = aicbic(Loglikehood,2,250);
res.aic(H,:)=a;
res.bic(H,:)=b;
clearvars a b Loglikehood Md1 sys R residual3 prediction3 

H=H+1;
end

% Autoregressive Integrated Moving Average (ARIMA)
for ar=5:5:25 %10, 15, 20, 25

ma=2; %2, 5
I=5; 
sys = arima(ar,I,ma);
[Md1 ,~,Loglikehood] = estimate(sys,dt_sa);
residual4 = infer(Md1,dt_tg);
prediction4 = dt_tg + residual4;

R = corrcoef(dt_tg,prediction4);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction4);
res.RMSE(H)=rms(dt_tg-prediction4);
res.model(H,:)=string(strcat('ARIMA(',num2str(ar),',',num2str(I),',',num2str(ma),')'));


% Goodness-of-Fit Checks
% Alaike or Bayesoan: to avoid the fitting is over-fit.
[a,b] = aicbic(Loglikehood,2,250);
res.aic(H,:)=a;
res.bic(H,:)=b;
clearvars a b Loglikehood Md1 sys R residual4 prediction4 

H=H+1;
end

for ar=5:5:25 %10, 15, 20, 25

ma=5; %2, 5
I=5; 
sys = arima(ar,I,ma);
[Md1 ,~,Loglikehood] = estimate(sys,dt_sa);
residual4 = infer(Md1,dt_tg);
prediction4 = dt_tg + residual4;

R = corrcoef(dt_tg,prediction4);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction4);
res.RMSE(H)=rms(dt_tg-prediction4);
res.model(H,:)=string(strcat('ARIMA(',num2str(ar),',',num2str(I),',',num2str(ma),')'));

% Goodness-of-Fit Checks
% Alaike or Bayesoan: to avoid the fitting is over-fit.
[a,b] = aicbic(Loglikehood,2,250);
res.aic(H,:)=a;
res.bic(H,:)=b;
clearvars a b Loglikehood Md1 sys R residual4 prediction4 

H=H+1;
end

for S=12:12:24
% Seasonal Autoregressive Integrated Moving-Average (SARIMA)
sys = arima('Constant',NaN,'ARLags',1:4,'D',0,'MALags',1:2,'SARLags',[12,24,36,48],'Seasonality',S,'SMALags',S,'Distribution','Gaussian');
[Md1 ,~,Loglikehood] = estimate(sys,dt_sa);
residual5 = infer(Md1,dt_tg);
prediction5 = dt_tg + residual5;


R = corrcoef(dt_tg,prediction5);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction5);
res.RMSE(H)=rms(dt_tg-prediction5);
res.model(H,:)=string(strcat('SARIMA(',num2str(S),')'));

% Goodness-of-Fit Checks
% Alaike or Bayesoan: to avoid the fitting is over-fit.
[a,b] = aicbic(Loglikehood,2,250);
res.aic(H,:)=a;
res.bic(H,:)=b;
clearvars a b Loglikehood Md1 sys R residual5 prediction5 

H=H+1;
end

% GARCH Model
GARCH_X1 = garch('Offset',0,'GARCHLags',1:14,'ARCHLags',50,'Distribution','Gaussian');
[GARCH_X1 ,~,Loglikehood]= estimate(GARCH_X1,dt_sa,'Display','off');
residual6 = infer(GARCH_X1,dt_tg);
prediction6 = dt_tg + residual6;

R = corrcoef(dt_tg,prediction6);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction6);
res.RMSE(H)=rms(dt_tg-prediction6);
res.model(H,:)="GARCH";

% Goodness-of-Fit Checks
% Alaike or Bayesoan: to avoid the fitting is over-fit.
[a,b] = aicbic(Loglikehood,2,250);
res.aic(H,:)=a;
res.bic(H,:)=b;
clearvars a b Loglikehood Md1 sys R residual6 prediction6 GARCH_X1
H=H+1;


%Glostan, Jagannathan and Runkle GARCH Model
% this give the same result as GRACH model in resampling
GJR_X1 = gjr('Offset',0,'GARCHLags',1:3,'ARCHLags',1,'LeverageLags',1,'Distribution','Gaussian');
[GJR_X1 ,~,Loglikehood] = estimate(GJR_X1,dt_sa,'Display','off');

residual7 = infer(GJR_X1,dt_tg);
prediction7 = dt_tg + residual7;

R = corrcoef(dt_tg,prediction7);
res.R(H) = R(1,2);
res.MAE(H)=mae(dt_tg-prediction7);
res.RMSE(H)=rms(dt_tg-prediction7);
res.model(H,:)="GJR";

% Goodness-of-Fit Checks
% Alaike or Bayesoan: to avoid the fitting is over-fit.
[a,b] = aicbic(Loglikehood,2,250);
res.aic(H,:)=a;
res.bic(H,:)=b;
clearvars a b Loglikehood GJR_X1 sys R residual7 prediction7 
H=H+1;

%% FIND THE BEST CASE
close all

for i=1:12
    
figure(i)
% colororder({'k','m'})

yyaxis left

plot(cell2mat(results{i}.AIC([1:17 32 33],:)),'^-r','MarkerFaceColor','r','DisplayName','AIC')
hold on
plot(cell2mat(results{i}.BIC([1:17 32 33],:)),'^-b','MarkerFaceColor','b','DisplayName','BIC')
set(gca,'YColor','m')


ylabel('criterion','FontSize',18,'FontWeight','bold');

yyaxis right

plot(cell2mat(results{i}.RMSE([1:17 32 33],:)),'o-k','MarkerFaceColor','k','DisplayName','RMSE')
plot(cell2mat(results{i}.MAE([1:17 32 33],:)),'o--k','MarkerFaceColor','k','DisplayName','MAE')
set(gca,'YColor','k')
xticks(1:1:19)
xlim([0 20])
ylabel('RMSE-MAE [cm]','FontSize',18,'FontWeight','bold');

xticklabels(results{i}.Model([1:17 32 33],:))

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);

legend show
ylim([4.5 10])

end

% subplot(2,2)

%The lowest AIC and BIC is best fit model
%% find group by catagory

for k=1:16
me=table();
[x,tm] = findgroups(sacut{k}.ym);

me.dtsa=splitapply(@mean, sacut{k}.sadt,x);
me.dttg=splitapply(@mean, sacut{k}.tgdt,x);
me.sadt_Linear=splitapply(@mean, sacut{k}.sadt_Linear,x);
me.t=tm;
sa_mean{k}=me;
clearvars me x tm

end

%Perform Differences
for k=1:16
sa_mean{k}.sadt_Differences=nan(height(sa_mean{k}),1);
sa_mean{k}.sadt_Differences(2:end)= diff(sa_mean{k}.sadt_Linear);

end

for k=1:16
% a=(sa_mean{k}.t);
sa_mean{k}.ye=round((sa_mean{k}.t)/100,0);
sa_mean{k}.mo=round((((sa_mean{k}.t)/100)-round((sa_mean{k}.t)/100,0))*100,0);
end

for k=1:16
sa_mean{k}.time=datetime(sa_mean{k}.ye,sa_mean{k}.mo,1);
end
%% plot acf pacf on monthly dt

close all

for k=1:14

figure(k)
subplot(3,3,1)
plot(sa_mean{k}.dtsa,'-','color',[0.4940 0.1840 0.5560])
hold on
p1 = polyfit(decyear(sa_mean{k}.t),sa_mean{k}.dtsa,1);
f1 = polyval(p1,decyear(sa_mean{k}.t));
tr1=fitlm(decyear(sa_mean{k}.t),sa_mean{k}.dtsa);
p(1)=plot(f1,'-k','LineWidth',2.5,'DisplayName',strcat('trend= ',num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'));
h=yline(mean(sa_mean{k}.dtsa),'--k',strcat('MDT=',num2str(mean(sa_mean{k}.dtsa),2),' [cm]'),'LineWidth',2);
h.LabelHorizontalAlignment='left';
h.FontSize=18;
legend(p(1),'FontSize',18,'Location','northwest')
legend boxoff
ylabel('DT [cm]')
xlabel([])
title('DT')

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);

subplot(3,3,2)
autocorr(sa_mean{k}.dtsa,'NumLags',20,'NumMA',1)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);
xlabel([])
ylabel('ACF')
title('DT')

% Partial ACF (PACF)
% Identify Series with Serial Correlation
% Determine whether an MA model is appropriate
% Identify significant AR lags for model identification.
subplot(3,3,3)
parcorr(sa_mean{k}.dtsa,'NumAR',1,'NumLags',20)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);
xlabel([])
ylabel('PACF')
title('DT')

subplot(3,3,4)
plot(sa_mean{k}.sadt_Linear,'-','color',[0.4940 0.1840 0.5560])
hold on
h=yline(mean(sa_mean{k}.sadt_Linear),'--k',strcat('Mean detrended=',num2str(mean(sa_mean{k}.sadt_Linear),2),' [cm]'),'LineWidth',2);
h.LabelHorizontalAlignment='left';
h.FontSize=18;
ylabel('DT [cm]')
xlabel([])
title('detrended DT')
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);

%linear
subplot(3,3,5)
autocorr(sa_mean{k}.sadt_Linear,'NumLags',20,'NumMA',1)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);
xlabel([])
ylabel('ACF')
title('detrended DT')

subplot(3,3,6)
parcorr(sa_mean{k}.sadt_Linear,'NumAR',1,'NumLags',20)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);
xlabel([])
ylabel('PACF')
title('detrended DT')


subplot(3,3,7)
plot(sa_mean{k}.sadt_Differences,'-','color',[0.4940 0.1840 0.5560])
hold on
h=yline(mean(sa_mean{k}.sadt_Differences,'omitnan'),'--k',strcat('Mean diff=',num2str(mean(sa_mean{k}.sadt_Differences,'omitnan'),2),' [cm]'),'LineWidth',2);
h.LabelHorizontalAlignment='left';
h.FontSize=18;
ylabel('DT [cm]')
xlabel('date')
title('diff detrended DT')
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);


%difference
subplot(3,3,8)
autocorr(sa_mean{k}.sadt_Differences,'NumLags',20,'NumMA',1)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);
ylabel('ACF')
title('diff detrended DT')

subplot(3,3,9)
parcorr(sa_mean{k}.sadt_Differences,'NumAR',1,'NumLags',20)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);
ylabel('PACF')
title('diff detrended DT')

end


%% SARIMA on monthgly mean data


% select TG data
clear res
clc

% select tg and associated data
k=4;

for i=1:10
AR{i}=[1:i];
end

for i=1:5
MA{i}=[1:i];
end


res=table();
H=1;
dt_sa=sa_mean{k}.sadt_Linear;
% SARIMA

    % Seasonal Autoregressive Integrated Moving-Average (SARIMA)

for S=12:6:24
    for i=1:10
        for m=1:5
            sys = arima('Constant',NaN,'ARLags',AR{i},'D',0,'MALags',MA{m},'SARLags',[6,12,24,36,48],'Seasonality',S,'SMALags',S,'Distribution','Gaussian');
            [Md1 ,~,Loglikehood] = estimate(sys,dt_sa);
            residual5 = infer(Md1,dt_tg);
            prediction5 = dt_tg + residual5;
            
            
            R = corrcoef(dt_tg,prediction5);
            res.R(H) = R(1,2);
            res.MAE(H)=mae(dt_tg-prediction5);
            res.RMSE(H)=rms(dt_tg-prediction5);
            res.model(H,:)=string(strcat('ARIMA(',num2str(max(AR{i})),'0',num2str(max(MA{m})),')','S(',num2str(S),')'));
            
            % Goodness-of-Fit Checks
            % Alaike or Bayesoan: to avoid the fitting is over-fit.
            [a,b] = aicbic(Loglikehood,2,250);
            res.aic(H,:)=a;
            res.bic(H,:)=b;
            clearvars a b Loglikehood Md1 sys R residual5 prediction5
            
            H=H+1;
        end
    end
end

%% prepare data

% TG data for forcasting comparision

%sort by date & monthly category
for k=1:16
tg{k} = sortrows(tg{k},'Time','ascend');
tg{k}.mo = month(tg{k}.Time);
tg{k}.ym=str2num([num2str(tg{k}.y),num2str(tg{k}.mo,'%0.2i')]);
end

for k=1:16
me=table();
[x,tm] = findgroups(tg{k}.ym);
me.dttg=splitapply(@mean, tg{k}.dt,x);
me.t=tm;
tg_mean{k}=me;
clearvars me x tm
end

for k=1:16
tg_mean{k}.ye=round((tg_mean{k}.t)/100,0);
tg_mean{k}.mo=round((((tg_mean{k}.t)/100)-round((tg_mean{k}.t)/100,0))*100,0);
tg_mean{k}.time=datetime(tg_mean{k}.ye,tg_mean{k}.mo,1);
tg_mean{k}.tg_detrend=detrend(tg_mean{k}.dttg);
end
%% Forcasting

%The lowest AIC and BIC is best fit model
k=11;
len=24;

t_sa=datetime(sa_mean{k}.ye,sa_mean{k}.mo,1);

for i=1:len
t_sa_forcast(i,1)=datetime(addtodate(datenum(t_sa(end)), i, 'month'),'ConvertFrom','datenum');
end

dt_sa=sa_mean{k}.sadt_Linear;

%change
dt_tg=detrend(tg{1, k}.dt,'omitnan');
date_tg=tg{1, k}.Time; 


S=12;
sys1 = arima('Constant',NaN,'ARLags',1:4,'D',0,'MALags',1:1,'SARLags',[12,24,36,48],'Seasonality',12,'SMALags',12,'Distribution','Gaussian');
Md1 = estimate(sys1,dt_sa);

sys2 = arima('Constant',NaN,'ARLags',1:4,'D',0,'MALags',1:1,'SARLags',[12,24,36,48],'Seasonality',6,'SMALags',6,'Distribution','Gaussian');
Md2 = estimate(sys2,dt_sa);

sys3 = arima('Constant',NaN,'ARLags',1:10,'D',0,'MALags',1:1,'SARLags',[12,24,36,48],'Seasonality',6,'SMALags',6,'Distribution','Gaussian');
Md3 = estimate(sys3,dt_sa);

forcast=table();
forcast.dt1 = forecast(Md1,len,'Y0',dt_sa);
forcast.dt2 = forecast(Md2,len,'Y0',dt_sa);
forcast.dt3 = forecast(Md3,len,'Y0',dt_sa);
forcast.t=t_sa_forcast;

for i=1:height(forcast)
forcast.tg_detrend(i)=tg_mean{k}.tg_detrend(tg_mean{k}.time==forcast.t(i));
forcast.tg(i)=tg_mean{k}.dttg(tg_mean{k}.time==forcast.t(i));
end


% retrend dt SA
% tr=fitlm(decyear(sacut{k}.t),sacut{k}.sadt);
% T=tr.Coefficients.Estimate(2);
% dt_Forecast1_retrend=retrend(dt_Forecast1,T);
% dt_Forecast2_retrend=retrend(dt_Forecast2,T);
% dt_Forecast3_retrend=retrend(dt_Forecast3,T);



    plot(tg_mean{k}.time,tg_mean{k}.tg_detrend,'-','color',[.7 .7 .7])
    hold on
    plot(t_sa,sa_mean{k}.sadt_Linear,'--*','color',[0.4940 0.1840 0.5560],'LineWidth',.5)

    %plot forcast
    h(1)=plot(forcast.t,forcast.dt1,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2,'DisplayName',strcat('S(12)ARIMA(4,0,1) ','RMSE=',num2str((rms(forcast.dt1-forcast.tg_detrend)),2),' [cm]'));
    h(2)=plot(forcast.t,forcast.dt2,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2,'DisplayName',strcat('S(6)ARIMA(4,0,1) ','RMSE=',num2str((rms(forcast.dt2-forcast.tg_detrend)),2),' [cm]'));
    h(3)=plot(forcast.t,forcast.dt3,'-','color',[0 0.4470 0.7410],'LineWidth',2,'DisplayName',strcat('S(6)ARIMA(10,0,1) ','RMSE=',num2str((rms(forcast.dt3-forcast.tg_detrend)),2),' [cm]'));
    
%     
%     plot(t_sa_forcast,dt_Forecast1_retrend,'-','color',[0.4660 0.6740 0.1880])
%     plot(t_sa_forcast,dt_Forecast2_retrend,'-','color',[0.6350 0.0780 0.1840])
%     plot(t_sa_forcast,dt_Forecast3_retrend,'-','color',[0 0.4470 0.7410])
    legend(h(1:3))
    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;
    ylabel('DT_d_e_t_e_r_e_n_d [cm]','FontSize',18,'FontWeight','bold');
    xlim([min(t_sa)-100,max(forcast.t)+100])
    pbaspect([1 .3 1])
 
    
    
