% gap filling
% load('G:\trend\data\alldata.mat')
load('F:\Analysis\BS\trend\data\alldata.mat')

for k=1:16
sa_fill{k}.sa=sa_fill{k}.sa_deseason+mean(sa_fill{k}.tg_deseason-sa_fill{k}.sa_deseason,'omitnan');
end

clearvars k 
%% Gap-filling

% determine the gaps size
gs=size(sa_fill{1}((191:218),:),1);


for k=1:16
    res{k}=table();
    gapsa{k}=table();

        %         forecasting
        dt_sa=sa_fill{k}.sa(5:190);
        % t_sa=sa_fill{k}.time(5:190);
        
        %     find the prediction length by acf
        %ACF
        [acf,malags,bounds] =  autocorr(dt_sa);
        % select lag of MA
        [n,~]=find(abs(acf)<=abs(bounds(1)),1);
        if malags(n-1)==0
            len=1;
        else
            len=malags(n-1);
        end
        
        H=1;
        HH=1;
        Hh=1;
        
        clearvars acf malags bounds
        
        % select Integrated lag  (I)
        if ~adftest(dt_sa)
            
            I=1;
            while ~adftest(diff(dt_sa,I))
                I=I+1;
            end
        else
            I=0;
        end
        
        
        
        while Hh<=ceil(gs/len)
            
            % select Integrated lag  (I)
            if ~adftest(dt_sa)
                
                I=1;
                while ~adftest(diff(dt_sa,I))
                    I=I+1;
                end
            else
                I=0;
            end
            
            
            
            %ACF
            [acf,malags,bounds] =  autocorr(dt_sa);
            % select lag of MA
            [n,~]=find(abs(acf)<=abs(bounds(1)),1);
            if malags(n-1)==0
                ma=1;
            else
                ma=malags(n-1);
            end
            
            %PACF
            [pacf,arlags,bounds] =  parcorr(dt_sa);
            % select lag of AR
            [n,~]=find(abs(pacf)<=abs(bounds(1)),1);
            if arlags(n-1)==0
                ar=1;
            else
                ar=arlags(n-1);
            end
            
            p=0+1;sys(p) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',12:6:48,'Seasonality',6,'SMALags',6,'Distribution','Gaussian');
            p=p+1;sys(p) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',12:6:48,'Seasonality',12,'SMALags',12,'Distribution','Gaussian');
            p=p+1;sys(p) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',12:6:48,'Seasonality',12,'SMALags',12:6:48,'Distribution','Gaussian');
            p=p+1;sys(p) = arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'SARLags',12:6:48,'Seasonality',6,'SMALags',12:6:48,'Distribution','Gaussian');
            p=p+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'Seasonality',6,'Distribution','Gaussian');
            p=p+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'Seasonality',12,'Distribution','Gaussian');
            p=p+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma);
            p=p+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar+1,'D',I,'MALags',1:ma+1);
            p=p+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar+2,'D',I,'MALags',1:ma+2);
            
            
            HH=1+(size(sys,2)*(Hh-1));
            
            name=table();
            
            for o=1:p
                name.a(o,1)=sys(o).Description;
            end
            
            res{k}.forename(HH:HH+size(sys,2)-1,1)=name.a;
            
            for j=1:p
                try
                    [~,~,Loglikehood] = estimate(sys(j),dt_sa,'Display','off');
                    [res{k}.foreaic(H,:),res{k}.forebic(H,:)] = aicbic(Loglikehood,2,250);
                    H=H+1;
                    clearvars Loglikehood
                catch ME
                    res{k}.foreaic(H,:)=NaN; res{k}.forebic(H,:)=NaN;
                    fprintf(strcat('FORECASTING did not work sys:',num2str(j),'for TG no:',num2str(k)));
                    H=H+1;
                    clearvars Loglikehood
                end
            end
            
            [n,~]=find(res{k}.foreaic(HH:HH+size(sys,2)-1,1)==min(res{k}.foreaic(HH:HH+size(sys,2)-1,1)));
            
            res{k}.foremodel(n+(size(sys,2)*(Hh-1)))=1;
            Md_signal = estimate(sys(n),dt_sa,'Display','off');
            
            pred(Hh,1)=forecast(Md_signal,len,'Y0',dt_sa);
            dt_sa=[dt_sa;pred(Hh)];
            clearvars Md_signal sys n ar ma
            Hh=Hh+1;
        end
        
        
        gapsa{k}.sadt_fore =pred;
        gapsa{k}.t=sa_fill{k}.time(191:218);
        gapsa{k}.tgdt=sa_fill{k}.tg_deseason(191:218);
        clearvars sys Md_signal ar ma I pacf arlags bounds acf malags detail j n a pred dt_sa H p Hh HH name dt_sa

        %         hindcasting
        dt_sa=flipud(sa_fill{k}.sa(219:258));
        % t_sa=sa_fill{k}.time(219:258);
        
        %Signal hindcast
        
        
        %     find the prediction length by acf
        %ACF
        [acf,malags,bounds] =  autocorr(dt_sa);
        % select lag of MA
        [n,~]=find(abs(acf)<=abs(bounds(1)),1);
        if malags(n-1)==0
            len=1;
        else
            len=malags(n-1);
        end
        clearvars acf malags bounds
        
        
        H=1;
        HH=1;
        Hh=1;
        
        while Hh<=ceil(gs/len)
            
            % select Integrated lag  (I)
            if ~adftest(dt_sa)
                
                I=1;
                while ~adftest(diff(dt_sa,I))
                    I=I+1;
                end
            else
                I=0;
            end
            
            
            %ACF
            [acf,malags,bounds] =  autocorr(dt_sa);
            % select lag of MA
            [n,~]=find(abs(acf)<=abs(bounds(1)),1);
            if malags(n-1)==0
                ma=1;
            else
                ma=malags(n-1);
            end
            
            %PACF
            [pacf,arlags,bounds] =  parcorr(dt_sa);
            % select lag of AR
            [n,~]=find(abs(pacf)<=abs(bounds(1)),1);
            if arlags(n-1)==0
                ar=1;
            else
                ar=arlags(n-1);
            end
            
            
            p=0+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'Seasonality',6,'Distribution','Gaussian');
            p=p+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma,'Seasonality',12,'Distribution','Gaussian');
            p=p+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma);
            p=p+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar+1,'D',I,'MALags',1:ma+1);
            p=p+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar+2,'D',I,'MALags',1:ma+2);
            p=p+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar,'D',I+1,'MALags',1:ma);
            p=p+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar,'D',I+2,'MALags',1:ma);
            p=p+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar+1,'D',I,'MALags',1:ma);
            p=p+1;sys(p)= arima('Constant',NaN,'ARLags',1:ar,'D',I,'MALags',1:ma+1);
            
            
            HH=1+(size(sys,2)*(Hh-1));
            name=table();
            
            for o=1:p
                name.a(o,1)=sys(o).Description;
            end
            
            res{k}.hindname(HH:HH+size(sys,2)-1,1)=name.a;

            
            for j=1:p
                try
                    [~,~,Loglikehood] = estimate(sys(j),dt_sa,'Display','off');
                    [res{k}.hindaic(H,:),res{k}.hindbic(H,:)] = aicbic(Loglikehood,2,250);
                catch ME
                    res{k}.hindaic(H,:)=NaN; res{k}.hindbic(H,:)=NaN;
                    fprintf(strcat('HINDCASTING did not work sys:',num2str(j),'for TG no:',num2str(k)));
                end
                H=H+1;
                clearvars Loglikehood
            end
            
            [n,~]=find(res{k}.hindaic(HH:HH+size(sys,2)-1,1)==min(res{k}.hindaic(HH:HH+size(sys,2)-1,1)));
            
            res{k}.hindmodel(n+(size(sys,2)*(Hh-1)))=1;
            Md_signal = estimate(sys(n),dt_sa,'Display','off');
            
            pred(Hh,1)=forecast(Md_signal,len,'Y0',dt_sa);
            dt_sa=[dt_sa;pred(Hh)];
            clearvars Md_signal sys n ar ma
            Hh=Hh+1;
        end
        
        % forcast table
        gapsa{k}.sadt_hind = flipud(pred);
        clearvars sys Md_signal ar ma I pacf arlags bounds acf malags detail j n a pred dt_sa H p Hh HH name o dt_sa
    
end

% clearvars res gapsa sys Md_signal ar ma I pacf arlags bounds acf malags detail j Hh H a pred dt_sa len HH o name prediction n p gs k ME

%% 

a=height(sa_fill{1}.sa(5:190));
b=height(sa_fill{1}.sa(219:258));

for k=1:16
    gapsa{k}.sadt_bi=nan(height(gapsa{k}.sadt_fore),1);
    for i=1:height(gapsa{k}.sadt_fore)
        gapsa{k}.sadt_bi(i,1)=((gapsa{k}.sadt_fore(i,1)*(a/(a+b)))+(gapsa{k}.sadt_hind(i,1)*(b/(a+b))))/2;
    end
    rmse_forcast(k,1)=rms(gapsa{k}.tgdt-gapsa{k}.sadt_fore,'omitnan');
    rmse_forcast(k,2)=rms(gapsa{k}.tgdt-gapsa{k}.sadt_hind,'omitnan');
    rmse_forcast(k,3)=rms(gapsa{k}.tgdt-gapsa{k}.sadt_bi,'omitnan');
    rmse_forcast(k,4)=rms(sa_fill{k}.dtsa-sa_fill{k}.dttg,'omitnan');

end


clearvars sys Md_signal ar ma I pacf arlags bounds acf malags detail j Hh H a pred dt_sa len HH o name prediction n p gs k ME


plot(rmse_forcast(:,1),'LineWidth',1.5,'DisplayName','Forecast','color',[0 0.4470 0.7410])
hold on
plot(rmse_forcast(:,2),'LineWidth',1.5,'DisplayName','Hindcast','color',[0.4940 0.1840 0.5560])
plot(rmse_forcast(:,3),'LineWidth',1.5,'DisplayName','Bi-Direction','color',[0.6350 0.0780 0.1840])
% plot(rmse_forcast(:,4),'LineWidth',1.5,'DisplayName','Original','color',[0.3010 0.7450 0.9330]	)
legend show
xticks([1:1:16])
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=27; ax.FontWeight='Bold'; grid on; ax.FontName='Times New Roman';
set(gca,'fontname','Times New Roman','FontSize',18);
xlim([0 17])

ylabel('RMSE [cm]','FontSize',20,'FontWeight','bold');
xlabel('TG ID','FontSize',20,'FontWeight','bold');


%% 
load('gapfill.mat')
H=1;  
model=table();
gs=size(sa_fill{1}((191:218),:),1);


for k=1:16
model.fore(H:H+gs-1,1)=res{k}.forename(res{k}.foremodel==1);
model.hind(H:H+gs-1,1)=res{k}.hindname(res{k}.hindmodel==1);
model.tgid(H:H+gs-1,1)=repmat(k,gs,1);
H=H+gs;
end

fore=unique(model.fore); %result: sys(3)
hind=unique(model.hind); % result: sys(5)

%% 



