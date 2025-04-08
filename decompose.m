%% decompositions
% EMD decomposition

% complete later
k=1;
dt_sa=sa_fill{k}.sa_deseason(5:190);
ti=sa_fill{k}.time(5:190);
dt_tg=sa_fill{k}.tg_deseason(5:190);



% EMD decompostion
[imf,residual,~] = emd(dt_sa,'MaxNumIMF',3);
[imf2,residual2,~] = emd(dt_tg,'MaxNumIMF',3);

forcast_imf{k}=table();
% forcast details

for i=1:size(imf,2)
    
    % select Integrated lag  (I)
    if ~adftest(imf(:,i))
        detail=imf(:,i);
        I=1;
        while ~adftest(diff(detail,I))
            I=I+1;
        end
    else
        I=0;
    end
    
          detail=imf(:,i);
  
    
    
    %ACF
    [acf,malags,bounds] =  autocorr(detail);
    % select lag of MA
    [n,~]=find(abs(acf)<=abs(bounds(1)),1);
    if malags(n-1)==0
        ma=1;
    else
        ma=malags(n-1);
    end
    
    %PACF
    [pacf,arlags,bounds]=parcorr(detail);
    % select lag of AR
    [n,~]=find(abs(pacf)<=abs(bounds(1)),1);
    if arlags(n-1)==0
        ar=1;
    else
        ar=malags(n-1);
    end
    
    
    if i~=3
        sys(1) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',[12,24,36,48],'Seasonality',6,'SMALags',6,'Distribution','Gaussian');
        sys(2) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',[12,24,36,48],'Seasonality',12,'SMALags',12,'Distribution','Gaussian');
        sys(3) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',[12,24,36,48],'Seasonality',12,'SMALags',[12,24,36,48],'Distribution','Gaussian');
        sys(4) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',[12,24,36,48],'Seasonality',6,'SMALags',[12,24,36,48],'Distribution','Gaussian');
        sys(5) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',12,'Seasonality',12,'SMALags',12,'Distribution','Gaussian');
        sys(6) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',6,'Seasonality',6,'SMALags',6,'Distribution','Gaussian');
        sys(7) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'Seasonality',6,'Distribution','Gaussian');
        sys(8) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'Seasonality',12,'Distribution','Gaussian');
        sys(9) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma);
        sys(10) = arima('Constant',NaN,'ARLags',1:ar+1,'D',I,'MALags',1:ma+1);
        sys(11) = arima('Constant',NaN,'ARLags',1:ar+1,'D',I,'MALags',1:ma);
        sys(12) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma+1);
        sys(13) = arima('Constant',1,'ARLags',1:ar,'D',I,'MALags',1:ma+1);
        sys(14) = arima('Constant',0,'ARLags',1:ar,'D',I,'MALags',1:ma+1);
    else
        sys(1) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'Seasonality',6,'Distribution','Gaussian');
        sys(2) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'Seasonality',12,'Distribution','Gaussian');
        sys(3) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma);
        sys(4) = arima('Constant',NaN,'ARLags',1:ar+1,'D',I,'MALags',1:ma+1);
        sys(5) = arima('Constant',NaN,'ARLags',1:ar+1,'D',I,'MALags',1:ma);
        sys(6) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma+1);
        sys(7) = arima('Constant',1,'ARLags',1:ar,'D',I,'MALags',1:ma+1);
        sys(8) = arima('Constant',0,'ARLags',1:ar,'D',I,'MALags',1:ma+1);
    end
    
    
    res=table();
    for j=1:size(sys,2)
        
        [~,~,Loglikehood] = estimate(sys(j),detail);
        [res.aic(j,:),res.bic(j,:)] = aicbic(Loglikehood,2,250);
        clearvars a b Loglikehood
    end
   
    [n,~]=find(res.aic==min(res.aic));
    Md_signal = estimate(sys(n),detail);
    
    len=28; %numel(sa_fill{1, 1}.sa_deseason(191:218))
    % forcast table
    forcast_imf{k}.t=sa_fill{k}.time(191:218);
    forcast_imf{k}.tgdt=sa_fill{k}.tg_deseason(191:218);
    forcast_imf{k}.("imf"+i) = forecast(Md_signal,len,'Y0',detail);
    
    clearvars sys Md_signal ar ma I pacf arlags bounds acf malags detail    
end



% forcast residuals
% select Integrated lag  (I)  
if ~adftest(residual)
    I=1;
    while ~adftest(diff(residual(1:i^2:end),I))
        I=I+1;
    detail=diff(detail,I);
    end
else
I=0; 
end

% select high information data only
detail=residual(1:2^i:end);


%ACF
[acf,malags,bounds] =  autocorr(detail);
% select lag of MA
[n,~]=find(abs(acf)<=abs(bounds(1)),1);
if malags(n-1)==0
    ma=1;
else
    ma=malags(n-1);
end

%PACF
[pacf,arlags,bounds] =  parcorr(detail);
% select lag of AR
[n,~]=find(abs(pacf)<=abs(bounds(1)),1);
if arlags(n-1)==0
    ar=1;
else
    ar=malags(n-1);
end

        sys(1) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'Seasonality',12,'Distribution','Gaussian');
        sys(2) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma);
        sys(3) = arima('Constant',NaN,'ARLags',1:ar+1,'D',I,'MALags',1:ma+1);
        sys(4) = arima('Constant',NaN,'ARLags',1:ar+1,'D',I,'MALags',1:ma);
        sys(5) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma+1);
        sys(6) = arima('Constant',1,'ARLags',1:ar,'D',I,'MALags',1:ma+1);
        sys(7) = arima('Constant',0,'ARLags',1:ar,'D',I,'MALags',1:ma+1);
   
    res=table();
    for j=1:size(sys,2)
        
        [~,~,Loglikehood] = estimate(sys(j),detail);
        [res.aic(j,:),res.bic(j,:)] = aicbic(Loglikehood,2,250);
        clearvars a b Loglikehood
    end
   
    [n,~]=find(res.aic==min(res.aic));
    Md_signal = estimate(sys(n),detail);

len=ceil(28/(2^i));
forcast=[detail;forecast(Md_signal,len,'Y0',detail)];

% forcastdiff=nan(numel(detail)+1+len,1);
% forcastdiff(1,1)=residual(1,1);
% forcastdiff(2:numel(detail)+1,1) = cumsum(detail)+residual(1,1);
% forcastdiff(numel(detail)+1:numel(detail)+1+len,1)=cumsum(forcast)+detail(1,1);

forcastint=interp1(1:2^i:length(forcast)*2^i,forcast,1:length(forcast)*2^i,'spline')';

forcast_imf{k}.res=forcastint(length(residual)+1:length(imf)+height(forcast_imf{k}));


for i=1:height(forcast_imf{k})
forcast_imf{k}.sadecom(i)=sum(forcast_imf{k}.imf1(i)+forcast_imf{k}.imf2(i)+forcast_imf{k}.imf3(i)-forcast_imf{k}.res(i));
end

[imf3,residual3,~] = emd(forcast_imf{k}.tgdt,'MaxNumIMF',3);
forcast_imf{k}.imftg1=imf3(:,1);
forcast_imf{k}.imftg2=imf3(:,2);
forcast_imf{k}.imftg3=imf3(:,3);
forcast_imf{k}.restg=residual3;

for i=1:height(forcast_imf{k})
    forcast_imf{1, 1}.dgdecom(i)=imf3(i,1)+imf3(i,2)+imf3(i,3);
end


for i=1:size(imf,2)+1
subplot(size(imf,2)+1,1,i)

if i~=4
    plot(1:186,imf(:,i))
    hold on
    plot(1:186,imf2(:,i))
    
    plot(187:186+28,forcast_imf{k}.("imftg"+i))
    plot(187:186+28,forcast_imf{k}.("imf"+i))

else
    plot(1:186,residual)
    hold on
    plot(1:186,residual2)    
    
    plot(187:186+28,forcast_imf{k}.restg)
    plot(187:186+28,forcast_imf{k}.res)
end

if i~=size(imf,2)+1
    xticklabels([])
end

xlim([0 216])

% legend show
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;



end



%% 
% EWT decomposition


k=1;
dt_sa=sa_fill{k}.sa_deseason(5:190);
ti=sa_fill{k}.time(5:190);
dt_tg=sa_fill{k}.tg_deseason(5:190);



[mra,cfs]  = ewt(dt_sa);
[mra2,cfs2]  = ewt(dt_tg);

forcast_mra{k}=table();
H=1;
resmra{k}=table();


for i=1:size(mra,2)
    
    % select Integrated lag  (I)
    if ~adftest(mra(:,i))
        I=1;
        while ~adftest(diff(detail,I))
            I=I+1;
        end
    else
        I=0;
    end
    

   detail=mra(:,i);
 
    
    %ACF
    [acf,malags,bounds] =  autocorr(detail);
    % select lag of MA
    [n,~]=find(abs(acf)<=abs(bounds(1)),1);
    if malags(n-1)==0
        ma=1;
    else
        ma=malags(n-1);
    end
    
    %PACF
    [pacf,arlags,bounds]=parcorr(detail);
    % select lag of AR
    [n,~]=find(abs(pacf)<=abs(bounds(1)),1);
    if arlags(n-1)==0
        ar=1;
    else
        ar=malags(n-1);
    end
    
    
   
        sys(1) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',[12,24,36,48],'Seasonality',6,'SMALags',6,'Distribution','Gaussian');
        sys(2) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',[12,24,36,48],'Seasonality',12,'SMALags',12,'Distribution','Gaussian');
        %sys(3) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',[12,24,36,48],'Seasonality',12,'SMALags',[12,24,36,48],'Distribution','Gaussian');
        sys(3) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',[12,24,36,48],'Seasonality',6,'SMALags',[12,24,36,48],'Distribution','Gaussian');
        sys(4) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',12,'Seasonality',12,'SMALags',12,'Distribution','Gaussian');
        sys(5) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',6,'Seasonality',6,'SMALags',6,'Distribution','Gaussian');
        sys(6) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'Seasonality',6,'Distribution','Gaussian');
        sys(7) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'Seasonality',12,'Distribution','Gaussian');
        sys(8) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma);
        sys(9) = arima('Constant',NaN,'ARLags',1:ar+1,'D',I,'MALags',1:ma+1);
        sys(10) = arima('Constant',NaN,'ARLags',1:ar+1,'D',I,'MALags',1:ma);
%         sys(11) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma+1);
%         sys(12) = arima('Constant',1,'ARLags',1:ar,'D',I,'MALags',1:ma+1);
%         sys(13) = arima('Constant',0,'ARLags',1:ar,'D',I,'MALags',1:ma+1);

    
    for j=1:size(sys,2)
        
        [~,~,Loglikehood] = estimate(sys(j),detail);
        [resmra{k}.aic(H,:),resmra{k}.bic(H,:)] = aicbic(Loglikehood,2,250);
        H=H+1;
        clearvars a b Loglikehood
    end
   
    Hh=i+size(sys,2)*(i-1);  
    [n,~]=find(resmra{k}.aic(Hh:H-1,1)==min(resmra{k}.aic(Hh:H-1,1)));
    Md_signal = estimate(sys(n),detail);
    
    len=28; %numel(sa_fill{1, 1}.sa_deseason(191:218))
    % forcast table
    forcast_mra{k}.t=sa_fill{k}.time(191:218);
    forcast_mra{k}.tgdt=sa_fill{k}.tg_deseason(191:218);
    forcast_mra{k}.("mra"+i) = forecast(Md_signal,len,'Y0',detail);
    
    clearvars sys Md_signal ar ma I pacf arlags bounds acf malags detail j
end
% clearvars resmra forcast_mra sys Md_signal ar ma I pacf arlags bounds acf malags detail j

for i=1:height(forcast_mra{k})
forcast_mra{k}.sadecom(i)=sum(forcast_mra{k}.mra1(i)+forcast_mra{k}.mra2(i)+forcast_mra{k}.mra3(i)-forcast_mra{k}.mra4(i));
end

%plot mra
for i=1:size(mra,2)+1
subplot(size(mra,2)+1,1,i)
if i==1
plot(ti,dt_sa,'--*','color',[0.4940 0.1840 0.5560],'LineWidth',.5)
hold on
plot(ti,dt_tg,'-','color',[.7 .7 .7])
hold off

xlim([min(ti)-100,max(ti)+100])
else
    plot(ti,mra(:,i-1))
    hold on
%    plot(ti,mra2(:,i-1))
    
end
    
if i==1
    title('DT_{deseasoned} [cm]')
    ylabel('DT','FontSize',18,'FontWeight','bold');
else
    ylabel(strcat('MRA',num2str(i)-1),'FontSize',18,'FontWeight','bold');
end

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;
xlim([min(ti)-1,max(ti)+28])

end
