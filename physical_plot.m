%% 
% tg_mean_m{1, 11}(1:2,:) = [];
% tg_mean_m{1, 10}(1:1296,:) = [];
% tg_mean_m{1, 9}(1:151,:) = [];
% tg_mean_m{1, 11}(697:end,:) = [];

%% 

close all

for i=7:13
[a,b]=xcorr(tg_mean_m{i}.dttg-mean(tg_mean_m{i}.dttg),'normalized');

plot(b,a,'DisplayName',strcat('TG_{id}:',num2str(i)),'LineWidth',1.5)
hold on
end
legend show
xlim([-150 150])
xlabel('Time Lag [month]')
ylabel('Autocorrelation of DT_T_G [cm]')
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';

%% decadal
close all
subplot(1,5,1)
scatter(ones(1,17),(1995:2011)',800,slide2_mean(1:17,1),"square",'filled')
hold on
scatter(repmat(2,1,17),(1995:2011)',800,slide2_mean(1:17,2),"square",'filled')
c=colorbar;
c.Label.String = 'mean sliding trend_{decadal} [mm/year]';

colormap(crameri('vik'))
xticks(1:2)
xticklabels(['SA';'TG'])
yticks(1995:2011);
yticklabels(deca)
% caxis([-15 15])
xlim([0 3])
ylim([1994 2012])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])


subplot(1,5,2)
scatter(ones(1,17),(1995:2011)',800,slide_ph1.th,"square",'filled')
c=colorbar;
xlabel('Temperature [\circC/yr]')

colormap(crameri('vik'))
yticks(1995:2011);
yticklabels([])
xticklabels([])

% caxis([-15 15])
xlim([0 2])
ylim([1994 2012])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])

subplot(1,5,3)
scatter(ones(1,17),(1995:2011)',800,slide_ph1.so,"square",'filled')
c=colorbar;
xlabel('Salinity [psu/yr]')

colormap(crameri('vik'))
yticks(1995:2011);
yticklabels([])
xticklabels([])

% caxis([-15 15])
xlim([0 2])
ylim([1994 2012])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])



subplot(1,5,4)
scatter(ones(1,17),(1995:2011)',800,slide_ph1.nao,"square",'filled')
c=colorbar;
xlabel('NAO/yr')

colormap(crameri('vik'))
yticks(1995:2011);
yticklabels([])
xticklabels([])

% caxis([-15 15])
xlim([0 2])
ylim([1994 2012])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])



subplot(1,5,5)
scatter(ones(1,17),(1995:2011)',800,slide_ph1.uo,"square",'filled')
hold on
scatter(repmat(2,1,17),(1995:2011)',800,slide_ph1.vo,"square",'filled')

c=colorbar;

colormap(crameri('vik'))
yticks(1995:2011);
yticklabels([])
xticks(1:2)
xticklabels(['EV';'NV'])

% caxis([-15 15])
xlim([0 3])
ylim([1994 2012])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])




%SA
r1=corrcoef(slide2_mean(1:17,1),slide_ph1.th);
r2=corrcoef(slide2_mean(1:17,1),slide_ph1.so);
r3=corrcoef(slide2_mean(1:17,1),slide_ph1.nao);
r4=corrcoef(slide2_mean(1:17,1),slide_ph1.uo);
r5=corrcoef(slide2_mean(1:17,1),slide_ph1.vo);


%TG
r11=corrcoef(slide2_mean(1:17,2),slide_ph1.th);
r12=corrcoef(slide2_mean(1:17,2),slide_ph1.so);
r13=corrcoef(slide2_mean(1:17,2),slide_ph1.nao);
r14=corrcoef(slide2_mean(1:17,2),slide_ph1.uo);
r15=corrcoef(slide2_mean(1:17,2),slide_ph1.vo);

r1(1,1)=r1(1,2)*100;
r1(2,1)=r2(1,2)*100;
r1(3,1)=r3(1,2)*100;
r1(4,1)=r4(1,2)*100;
r1(5,1)=r5(1,2)*100;


r1(1,2)=r11(1,2)*100;
r1(2,2)=r12(1,2)*100;
r1(3,2)=r13(1,2)*100;
r1(4,2)=r14(1,2)*100;
r1(5,2)=r15(1,2)*100;

r1

rms(slide2_mean(:,2)-slide2_mean(:,1))


%% interdecadal (y=5)
% close all
figure(2)
subplot(1,5,1)
scatter(ones(1,22),(1995:2016)',800,slide3_mean(1:22,1),"square",'filled')
hold on
scatter(repmat(2,1,22),(1995:2016)',800,slide3_mean(1:22,2),"square",'filled')
c=colorbar;
c.Label.String = 'mean sliding trend_{inter-decadal} [mm/year]';

colormap(crameri('vik'))
xticks(1:2)
xticklabels(['SA';'TG'])
yticks(1995:2016);
yticklabels(interd)
% caxis([-45 70])
xlim([0 3])
ylim([1994 2017])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])




subplot(1,5,2)
scatter(ones(1,22),(1995:2016)',800,slide_ph2.th,"square",'filled')
c=colorbar;
xlabel('Temperature [\circC/yr]')
colormap(crameri('vik'))
yticks(1995:2016);
ylim([1994 2017])

yticklabels([])
xticklabels([])

% caxis([-15 15])
xlim([0 2])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])

subplot(1,5,3)
scatter(ones(1,22),(1995:2016)',800,slide_ph2.so,"square",'filled')
c=colorbar;
xlabel('Salinity [psu/yr]')

colormap(crameri('vik'))
yticks(1995:2016);
ylim([1994 2017])

yticklabels([])
xticklabels([])

% caxis([-15 15])
xlim([0 2])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])



subplot(1,5,4)
scatter(ones(1,22),(1995:2016)',800,slide_ph2.nao,"square",'filled')
c=colorbar;
xlabel('NAO/yr')

colormap(crameri('vik'))
yticks(1995:2016);
ylim([1994 2017])

yticklabels([])
xticklabels([])

% caxis([-15 15])
xlim([0 2])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])




subplot(1,5,5)
scatter(ones(1,22),(1995:2016)',800,slide_ph2.uo,"square",'filled')
hold on
scatter(repmat(2,1,22),(1995:2016)',800,slide_ph2.vo,"square",'filled')

c=colorbar;

colormap(crameri('vik'))
yticks(1995:2016);
yticklabels([])
xticks(1:3)
xticklabels(['EV';'NV'])

% caxis([-15 15])
xlim([0 3])
ylim([1994 2017])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])




%SA
r1=corrcoef(slide3_mean(1:22,1),slide_ph2.th);
r2=corrcoef(slide3_mean(1:22,1),slide_ph2.so);
r3=corrcoef(slide3_mean(1:22,1),slide_ph2.nao);
r4=corrcoef(slide3_mean(1:22,1),slide_ph2.uo);
r5=corrcoef(slide3_mean(1:22,1),slide_ph2.vo);

%TG
r11=corrcoef(slide3_mean(1:22,2),slide_ph2.th);
r12=corrcoef(slide3_mean(1:22,2),slide_ph2.so);
r13=corrcoef(slide3_mean(1:22,2),slide_ph2.nao);
r14=corrcoef(slide3_mean(1:22,2),slide_ph2.uo);
r15=corrcoef(slide3_mean(1:22,2),slide_ph2.vo);

r2(1,1)=r1(1,2)*100;
r2(2,1)=r2(1,2)*100;
r2(3,1)=r3(1,2)*100;
r2(4,1)=r4(1,2)*100;
r2(5,1)=r5(1,2)*100;


r2(1,2)=r11(1,2)*100;
r2(2,2)=r12(1,2)*100;
r2(3,2)=r13(1,2)*100;
r2(4,2)=r14(1,2)*100;
r2(5,2)=r15(1,2)*100;


r2

rms(slide3_mean(:,2)-slide3_mean(:,1))

%% interdecadal (y=7)
% close all
figure(3)
subplot(1,5,1)
scatter(ones(1,20),(1995:2014)',800,slide4_mean(1:20,1),"square",'filled')
hold on
scatter(repmat(2,1,20),(1995:2014)',800,slide4_mean(1:20,2),"square",'filled')
c=colorbar;
c.Label.String = 'mean sliding trend_{inter-decadal} [mm/year]';

colormap(crameri('vik'))
xticks(1:2)
xticklabels(['SA';'TG'])
yticks(1995:2014);
yticklabels(interd2)
% caxis([-45 70])
xlim([0 3])
ylim([1994 2015])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])




subplot(1,5,2)
scatter(ones(1,20),(1995:2014)',800,slide_ph3.th,"square",'filled')
c=colorbar;
xlabel('Temperature [\circC/yr]')
colormap(crameri('vik'))
yticks(1995:2014);
ylim([1994 2015])

yticklabels([])
xticklabels([])

% caxis([-15 15])
xlim([0 2])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])

subplot(1,5,3)
scatter(ones(1,20),(1995:2014)',800,slide_ph3.so,"square",'filled')
c=colorbar;
xlabel('Salinity [psu/yr]')

colormap(crameri('vik'))
yticks(1995:2014);
ylim([1994 2015])

yticklabels([])
xticklabels([])

% caxis([-15 15])
xlim([0 2])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])



subplot(1,5,4)
scatter(ones(1,20),(1995:2014)',800,slide_ph3.nao,"square",'filled')
c=colorbar;
xlabel('NAO/yr')

colormap(crameri('vik'))
yticks(1995:2014);
ylim([1994 2015])

yticklabels([])
xticklabels([])

% caxis([-15 15])
xlim([0 2])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])





subplot(1,5,5)
scatter(ones(1,20),(1995:2014)',800,slide_ph3.uo,"square",'filled')
hold on
scatter(repmat(2,1,20),(1995:2014)',800,slide_ph3.vo,"square",'filled')

c=colorbar;

colormap(crameri('vik'))
yticks(1995:2014);
yticklabels([])
xticks(1:2)
xticklabels(['EV';'NV'])

% caxis([-15 15])
xlim([0 3])
ylim([1994 2015])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])





%SA
r1=corrcoef(slide4_mean(1:20,1),slide_ph3.th);
r2=corrcoef(slide4_mean(1:20,1),slide_ph3.so);
r3=corrcoef(slide4_mean(1:20,1),slide_ph3.nao);
r4=corrcoef(slide4_mean(1:20,1),slide_ph3.uo);
r5=corrcoef(slide4_mean(1:20,1),slide_ph3.vo);

%TG
r11=corrcoef(slide4_mean(1:20,2),slide_ph3.th);
r12=corrcoef(slide4_mean(1:20,2),slide_ph3.so);
r13=corrcoef(slide4_mean(1:20,2),slide_ph3.nao);
r14=corrcoef(slide4_mean(1:20,2),slide_ph3.uo);
r15=corrcoef(slide4_mean(1:20,2),slide_ph3.vo);


r3(1,1)=r1(1,2)*100;
r3(2,1)=r2(1,2)*100;
r3(3,1)=r3(1,2)*100;
r3(4,1)=r4(1,2)*100;
r3(5,1)=r5(1,2)*100;

r3(1,2)=r11(1,2)*100;
r3(2,2)=r12(1,2)*100;
r3(3,2)=r13(1,2)*100;
r3(4,2)=r14(1,2)*100;
r3(5,2)=r15(1,2)*100;


r3
rms(slide4_mean(:,2)-slide4_mean(:,1))



%% SLA

mdt_all=table();
for i=1:13
    k=1;
for y=1995:2022
mdt_all.sa(k,1)=mean(mean_all{i}.sa(year(mean_all{i}.time)==y));
mdt_all.tg(k,1)=mean(mean_all{i}.tg(year(mean_all{i}.time)==y));
mdt_all.y(k,1)=y;
k=k+1;
end
end

for i=1:27
ph.sla_sa(i,1)=ph.sa(i)-mean(ph.sa);
ph.sla_tg(i,1)=ph.tg(i)-mean(ph.sa);
end

%% plot NAO


r1=corrcoef(ph.sa,ph.nao);
rr1=r1(1,2)*100;

r2=corrcoef(ph.tg,ph.nao);
rr2=r2(1,2)*100;


plot(ph.y,ph.nao,'LineWidth',1.5,'color',[0.4660 0.6740 0.1880],'DisplayName','NAO')
ylabel('NAO')

yyaxis right
plot(ph.y,ph.sla_sa,'LineWidth',1.5,'color',[0.4940 0.1840 0.5560],'DisplayName','SA_{SLA}');
hold on
plot(ph.y,ph.sla_tg,'LineWidth',1.5,'color',[0.4940 0.1840 0.5560],'DisplayName','SA_{TG}');
ylabel('SLA [cm]')

legend show
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
grid minor
grid on
xlim([1990 2025])
pbaspect([1 .2 .5])