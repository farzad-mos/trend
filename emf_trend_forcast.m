%forcasting by EMD 
dt_sa=sa_mean{k}.sadt_Linear;
dt_tg=tg_mean{k}.tg_detrend;

[imf,residual,info] = emd(dt_sa);
% [imf,residual,info] = emd(dt_sa,'MaxNumIMF',7);

% [hs,f,t,imfinsf,imfinse] = hht(imf);

[mra,cfs]  = ewt(dt_sa);
% hht(mra)

figure(1)
for i=1:size(imf,2)+2
subplot(size(imf,2)+2,1,i)

if i==1
plot(sa_mean{k}.time,dt_sa,'--*','color',[0.4940 0.1840 0.5560],'LineWidth',.5)
hold on
plot(tg_mean{k}.time,tg_mean{k}.tg_detrend,'-','color',[.7 .7 .7])
hold off

xlim([min(sa_mean{k}.time)-100,max(sa_mean{k}.time)+100])

elseif i==size(imf,2)+2
    plot(sa_mean{k}.time,residual)
else
    plot(sa_mean{k}.time,imf(:,i-1))
end

if i==1
    title('DT_d_e_t_e_r_e_n_d [cm]')
    ylabel('DT','FontSize',18,'FontWeight','bold');
elseif i==size(imf,2)+2
    ylabel('Residuals','FontSize',18,'FontWeight','bold');
else
    ylabel(strcat('IMF',num2str(i)-1),'FontSize',18,'FontWeight','bold');
end

if i~=size(imf,2)+2
    xticklabels([])
end


ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;
xlim([min(sa_mean{k}.time)-100,max(sa_mean{k}.time)+100])

end


figure(2)
for i=1:size(mra,2)+1
subplot(size(mra,2)+1,1,i)
if i==1
plot(sa_mean{k}.time,dt_sa,'--*','color',[0.4940 0.1840 0.5560],'LineWidth',.5)
hold on
plot(tg_mean{k}.time,tg_mean{k}.tg_detrend,'-','color',[.7 .7 .7])
hold off

xlim([min(sa_mean{k}.time)-100,max(sa_mean{k}.time)+100])
else
    plot(sa_mean{k}.time,mra(:,i-1))
end
    
if i==1
    title('DT_d_e_t_e_r_e_n_d [cm]')
    ylabel('DT','FontSize',18,'FontWeight','bold');
else
    ylabel(strcat('MRA',num2str(i)-1),'FontSize',18,'FontWeight','bold');
end

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;
xlim([min(sa_mean{k}.time)-100,max(sa_mean{k}.time)+100])

end


%plot ACF PACF for each IMF

%Perform Differences for IMF1
imf1_dtsa= diff(imf(:,1));


for i=1:size(imf,2)

figure(i)
subplot(3,1,1)
plot(sa_mean{k}.time,imf(:,i),'-')
hold on
p1 = polyfit(decyear(sa_mean{k}.time),imf(:,i),1);
f1 = polyval(p1,decyear(sa_mean{k}.time));
tr1=fitlm(decyear(sa_mean{k}.time),imf(:,i));
p(1)=plot(sa_mean{k}.time,f1,'-k','LineWidth',2.5,'DisplayName',strcat('trend= ',num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'));
h=yline(mean(imf(:,i)),'--k',strcat('Mean=',num2str(mean(imf(:,i)),2),' [cm]'),'LineWidth',2);
h.LabelHorizontalAlignment='left';
h.FontSize=18;
legend(p(1),'FontSize',18,'Location','northwest')
legend boxoff
ylabel(strcat('IMF',num2str(i)))
xlabel('Date')

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; 
set(gca,'fontname','Times New Roman','FontSize',18);

subplot(3,1,2)
autocorr(imf(:,i),'NumLags',20,'NumMA',1)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; 
set(gca,'fontname','Times New Roman','FontSize',18);
% xlabel([])
ylabel('ACF')
title([])

% Partial ACF (PACF)
subplot(3,1,3)
parcorr(imf(:,i),'NumAR',1,'NumLags',20)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; 
set(gca,'fontname','Times New Roman','FontSize',18);
% xlabel([])
ylabel('PACF')
title([])
end

%% ACF PACF of diff IMF1

%Perform Differences for IMF1
imf1_dtsa= diff(imf(:,1));

%plot ACF PACF for each IMF
figure(i)
subplot(3,1,1)
plot(sa_mean{k}.time(2:end),imf1_dtsa,'-')
hold on
p1 = polyfit(decyear(sa_mean{k}.time(2:end)),imf1_dtsa,1);
f1 = polyval(p1,decyear(sa_mean{k}.time(2:end)));
tr1=fitlm(decyear(sa_mean{k}.time(2:end)),imf1_dtsa);
p(1)=plot(sa_mean{k}.time(2:end),f1,'-k','LineWidth',2.5,'DisplayName',strcat('trend= ',num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'));
h=yline(mean(imf1_dtsa),'--k',strcat('Mean=',num2str(mean(imf1_dtsa),2),' [cm]'),'LineWidth',2);
h.LabelHorizontalAlignment='left';
h.FontSize=18;
legend(p(1),'FontSize',18,'Location','northwest')
legend boxoff
ylabel(strcat('IMF',num2str(i),'_d_i_f_f'))
xlabel('Date')

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; 
set(gca,'fontname','Times New Roman','FontSize',18);

subplot(3,1,2)
autocorr(imf1_dtsa,'NumLags',20,'NumMA',1)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; 
set(gca,'fontname','Times New Roman','FontSize',18);
% xlabel([])
ylabel('ACF')
title([])

% Partial ACF (PACF)
subplot(3,1,3)
parcorr(imf1_dtsa,'NumAR',1,'NumLags',20)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; 
set(gca,'fontname','Times New Roman','FontSize',18);
% xlabel([])
ylabel('PACF')
title([])



%% ACF PACF residuals

subplot(3,1,1)
plot(sa_mean{k}.time,residual,'-')
hold on
p1 = polyfit(decyear(sa_mean{k}.time),residual,1);
f1 = polyval(p1,decyear(sa_mean{k}.time));
tr1=fitlm(decyear(sa_mean{k}.time),residual);
p(1)=plot(sa_mean{k}.time,f1,'-k','LineWidth',2.5,'DisplayName',strcat('trend= ',num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'));
h=yline(mean(residual),'--k',strcat('Mean=',num2str(mean(residual),2),' [cm]'),'LineWidth',2);
h.LabelHorizontalAlignment='left';
h.FontSize=18;
legend(p(1),'FontSize',18,'Location','northwest')
legend boxoff
ylabel(strcat('Residuals'))
xlabel('Date')

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; 
set(gca,'fontname','Times New Roman','FontSize',18);

subplot(3,1,2)
autocorr(residual,'NumLags',20,'NumMA',1)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; 
set(gca,'fontname','Times New Roman','FontSize',18);
% xlabel([])
ylabel('ACF')
title([])

% Partial ACF (PACF)
subplot(3,1,3)
parcorr(residual,'NumAR',1,'NumLags',20)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; 
set(gca,'fontname','Times New Roman','FontSize',18);
% xlabel([])
ylabel('PACF')
title([])


%% EMD forcasting

S=12;

% signal Model
sys=arima('Constant',NaN,'ARLags',1:4,'D',0,'MALags',1:1,'SARLags',[12,24,36,48],'Seasonality',12,'SMALags',12,'Distribution','Gaussian');

% IMFs Model
% sys1 = arima('Constant',NaN,'ARLags',1:7,'D',0,'MALags',1:2,'SARLags',[12,24,36,48],'Seasonality',S,'SMALags',S,'Distribution','Gaussian'); % diff
sys1 = arima('Constant',NaN,'ARLags',1:4,'D',0,'MALags',1:2,'SARLags',[12,24,36,48],'Seasonality',S,'SMALags',S,'Distribution','Gaussian');
sys2 = arima('Constant',NaN,'ARLags',1:8,'D',0,'MALags',1:2,'SARLags',[12,24,36,48],'Seasonality',S,'SMALags',S,'Distribution','Gaussian');
sys3 = arima('Constant',NaN,'ARLags',1:5,'D',0,'MALags',1:4,'SARLags',[12,24,36,48],'Seasonality',S,'SMALags',S,'Distribution','Gaussian');
sys4 = arima(6,0,0);
sys5 = arima(7,0,0);
sys6 = arima(5,0,0);

% Residuals Model
sys_r = arima(8,0,0);

% Signal forcasting
Md_signal= estimate(sys,dt_sa);

% EMD forcasting
% Md_imf1= estimate(sys1,imf1_dtsa); % diff has been used instead of imf(1)
Md_imf1= estimate(sys1,imf(:,1));

Md_imf2= estimate(sys2,imf(:,2));
Md_imf3= estimate(sys3,imf(:,3));
Md_imf4= estimate(sys4,imf(:,4));
Md_imf5= estimate(sys5,imf(:,5));
Md_imf6= estimate(sys6,imf(:,6));

% Residuals Forcasting
Md_residuals= estimate(sys_r,residual);


% forcast table
len=24;

t_sa=datetime(sa_mean{k}.ye,sa_mean{k}.mo,1);

for i=1:len
t_sa_forcast(i,1)=datetime(addtodate(datenum(t_sa(end)), i, 'month'),'ConvertFrom','datenum');
end

forcast_EMD=table();
forcast_EMD.t=t_sa_forcast;

%Signal forcast
forcast_EMD.dt = forecast(Md_signal,len,'Y0',dt_sa);


%EMD

%IMFs forcast
% diff
%  forcastdiff=forecast(Md_imf1,len-1,'Y0',imf1_dtsa);
%  forcast_EMD.dt1(1)=imf(1,1);
%  forcast_EMD.dt1(2:end) = cumsum(forcastdiff)+imf(1,1);

 forcast_EMD.dt11=forecast(Md_imf1,len,'Y0',imf1_dtsa);

forcast_EMD.dt2 = forecast(Md_imf2,len,'Y0',imf(:,2));
forcast_EMD.dt3 = forecast(Md_imf3,len,'Y0',imf(:,3));
forcast_EMD.dt4 = forecast(Md_imf4,len,'Y0',imf(:,4));
forcast_EMD.dt5 = forecast(Md_imf5,len,'Y0',imf(:,5));
forcast_EMD.dt6 = forecast(Md_imf6,len,'Y0',imf(:,6));

% residuals forcast
forcast_EMD.residuals = forecast(Md_residuals,len,'Y0',residual);


% forcast_EMD.dtEMD=forcast_EMD.dt1+forcast_EMD.dt2+forcast_EMD.dt3+forcast_EMD.dt4+forcast_EMD.dt5+forcast_EMD.dt6;
forcast_EMD.dtEMD2=forcast_EMD.dt11+forcast_EMD.dt2+forcast_EMD.dt3+forcast_EMD.dt4+forcast_EMD.dt5+forcast_EMD.dt6;
forcast_EMD.dtEMD3=forcast_EMD.dt11+forcast_EMD.dt2+forcast_EMD.dt3+forcast_EMD.dt4+forcast_EMD.dt5+forcast_EMD.dt6+forcast_EMD.residuals;


% add TG data at forcasting time
for i=1:height(forcast_EMD)
forcast_EMD.tg_detrend(i)=tg_mean{k}.tg_detrend(tg_mean{k}.time==forcast_EMD.t(i));
forcast_EMD.tg(i)=tg_mean{k}.dttg(tg_mean{k}.time==forcast_EMD.t(i));
end


h(1)=plot(forcast_EMD.t,forcast_EMD.dt,'-b','LineWidth',2,'DisplayName',strcat('signal_f_o_r_c_a_s_t ','RMSE=',num2str((rms(forcast_EMD.dt-forcast_EMD.tg_detrend)),3),' [cm]'));
hold on
% h(2)=plot(forcast_EMD.t,forcast_EMD.dtEMD2,'-r','LineWidth',2,'DisplayName',strcat('EMD_f_o_r_c_a_s_t ','RMSE=',num2str((rms(forcast_EMD.dtEMD2-forcast_EMD.tg_detrend)),4),' [cm]'));
% h(3)=plot(forcast_EMD.t,forcast_EMD.dtEMD,'-g','LineWidth',2,'DisplayName',strcat('EMD forcast','RMSE=',num2str((rms(forcast_EMD.dtEMD-forcast_EMD.tg_detrend)),2),' [cm]'));
h(4)=plot(forcast_EMD.t,forcast_EMD.dtEMD3,'-g','LineWidth',2,'DisplayName',strcat('EMD_f_o_r_c_a_s_t ','RMSE=',num2str((rms(forcast_EMD.dtEMD3-forcast_EMD.tg_detrend)),4),' [cm]'));
h(2)=plot(forcast_EMD.t,b(187:210),'-r','LineWidth',2,'DisplayName',strcat('EMD_f_o_r_c_a_s_t ','RMSE=',num2str((rms(b(187:210)-forcast_EMD.tg_detrend)),4),' [cm]'));

plot(t_sa,sa_mean{k}.sadt_Linear,'--*','color',[0.4940 0.1840 0.5560],'LineWidth',.5,'DisplayName','SA')

plot(tg_mean{k}.time,tg_mean{k}.tg_detrend,'-','color',[.7 .7 .7],'DisplayName','TG')

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;
ylabel('DT_d_e_t_e_r_e_n_d [cm]','FontSize',18,'FontWeight','bold');
xlim([min(t_sa)-100,max(forcast_EMD.t)+100])
pbaspect([1 .3 1])
% legend(h(1:2))
legend show


%% EWT forcasting

%
Md_mra1= estimate(sys1,mra(:,1));
Md_mra2= estimate(sys1,mra(:,2));
Md_mra3= estimate(sys1,mra(:,3));

%% 

figure(1)
for i=1:5
subplot(5,1,i)

if i==1
plot(sa_mean{k}.time,dt_sa,'--*','color',[0.4940 0.1840 0.5560],'LineWidth',.5)
hold on
plot(tg_mean{k}.time,tg_mean{k}.tg_detrend,'-','color',[.7 .7 .7])
hold off

xlim([min(sa_mean{k}.time)-100,max(sa_mean{k}.time)+100])

else
    plot(sa_mean{k}.time,dat(:,i-1))
end

if i==1
    ylabel('DT_d_e_t_e_r_e_n_d [cm]','FontSize',18,'FontWeight','bold');

else
    ylabel(strcat('D'),'FontSize',18,'FontWeight','bold');
end

if i~=5
    xticklabels([])
end


ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;
xlim([min(sa_mean{k}.time)-100,max(sa_mean{k}.time)+100])

end

