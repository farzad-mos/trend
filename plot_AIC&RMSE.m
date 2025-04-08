
subplot(2,1,2)

yyaxis left
h(1)=plot(res.aic,'^-r','MarkerFaceColor','r','DisplayName','AIC');
hold on
h(2)=plot(res.bic,'^-b','MarkerFaceColor','b','DisplayName','BIC');
yline(mean(res.aic),'-m','LineWidth',1.5)
yline(mean(res.bic),'--m','LineWidth',1.5)
set(gca,'YColor','m')
set(gca, 'YGrid', 'off', 'XGrid', 'off')

yyaxis right
h(3)=plot(res.RMSE,'o-k','MarkerFaceColor','k','DisplayName','RMSE');
hold on
h(4)=plot(res.MAE,'o--k','MarkerFaceColor','k','DisplayName','MAE');
xticks(1:1:71)
yline(mean(res.RMSE),'-k','LineWidth',1.5)
yline(mean(res.MAE),'--k','LineWidth',1.5)

grid on;
xlim([0 9])
xticklabels(res.model)
set(gca,'YColor','k')

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; ax.YGrid = 'on';  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',18);

legend(h(1:4))
%% 

S=6;
% i=4;
m=1;

res=table();
H=1;
dt_sa=sa_mean{k}.sadt_Linear;


    for i=1:9
            sys = arima('Constant',NaN,'ARLags',AR{i},'D',0,'MALags',MA{m},'SARLags',[6,12,24,36,48],'Seasonality',S,'SMALags',S,'Distribution','Gaussian');
            [Md1 ,~,Loglikehood] = estimate(sys,dt_sa);
            residual5 = infer(Md1,dt_tg);
            prediction5 = dt_tg + residual5;
            
            
            R = corrcoef(dt_tg,prediction5);
            res.R(H) = R(1,2);
            res.MAE(H)=mae(dt_tg-prediction5);
            res.RMSE(H)=rms(dt_tg-prediction5);
            res.model(H,:)=string(strcat('ARIMA(',num2str(max(AR{i})),',0,',num2str(max(MA{m})),')','S(',num2str(S),')'));
            
            % Goodness-of-Fit Checks
            % Alaike or Bayesoan: to avoid the fitting is over-fit.
            [a,b] = aicbic(Loglikehood,2,250);
            res.aic(H,:)=a;
            res.bic(H,:)=b;
            clearvars a b Loglikehood Md1 sys R residual5 prediction5
            
            H=H+1;
    end
