load('F:\Analysis\BS\trend\data\sa\all\SA_1.mat')

%% remove gross errors & outliers

clearvars  h2 a1 a2 e1 e2 h1


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

clearvars i k s out m j id


%% find within 50km distance data at TG
% extract points within 50 km distance of TG
%for said=6 (JA3) the radius is 70 km

for i=1:height(TGinfo)
    [lat,lon] = scircle1(TGinfo.Lat(i),TGinfo.Lon(i),km2deg(60),[],[],[],20);
    in = inpolygon(sa.lat,sa.lon,lat,lon);
    sa_tg{i}=sa(in==1,:);
    % [~,d] = dsearchn([Lat Lon],[sa.lat sa.lon]);
    clearvars i in lat lon
end

clear i

%% TG treatment
% combine all TG tables
tg={TG2,TG3,TG4,TG5,TG6,TG7,TG8,TG9,TG10,TG11,TG12,TG13,TG14};
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
    tg{i}.dt2=tg{i}.LU+tg{i}.zerotomean;

    clearvars T ia
end

i=9; %DK tg treatment
nontomean=(0.29541*0.7*(sind(TGinfo.Lat(i)).^2-sind(NAPlat).^2)+0.00042*(sind(TGinfo.Lat(i)).^4-sind(NAPlat).^4))*100; %add to TG data
tg{i}.dt2=tg{i}.LU+nontomean;

clear NAPlat i 

%% Plot VLM and geoid rise
% figure 2
load('Z:\data\raw_data\geoid_data\nkg2016lu-with-readme\NKG2016LU_abs.dat')
load('Z:\data\raw_data\geoid_data\nkg2016lu-with-readme\NKG2016LU_lev.dat')

xx = linspace(min(NKG2016LUabs.lat),max(NKG2016LUabs.lat),20);
yy = linspace(min(NKG2016LUabs.lon),max(NKG2016LUabs.lon),20);
[X,Y] = meshgrid(xx,yy);

Z = griddata(NKG2016LUabs.lat,NKG2016LUabs.lon,NKG2016LUabs.lu,X,Y);
Z2 = griddata(NKG2016LUabs.lat,NKG2016LUabs.lon,NKG2016LUabs.lu-NKG2016LUlev.lu,X,Y); %geoid rise
% plot VLM and GC
close all

contourfm(X,Y,Z,200,'edgecolor','none')
hold on
[M,C] = contour(Y,X,Z2,'-r','ShowText','on');
C.LineWidth = 1;
clabel(M,C,'FontSize',8,'Color','red')


plot(bs(:,2),bs(:,1),'.','color',[.5 .5 .5])
c=colorbar;
c.Label.String = 'VLM [mm/year]';
caxis([-4 10])
% colormap((turbo))
colormap(crameri('roma'))

for i=1:length(tg)
    plot(TGinfo.Lon(i),TGinfo.Lat(i),'Marker', '^','Color','k','MarkerFaceColor','k','MarkerSize',20)
end

text(TGinfo.Lon,TGinfo.Lat-0.3, TGinfo.TG,'FontSize',17,'Color','k','FontWeight','Bold');
xlim([9 30])
ylim([53 66.2])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; ax.FontName='Times New Roman';
c.FontSize=24;

%% remove geoid rise effect
epoch = 2000;

for i=1:length(tg)

    TGinfo.GC(i) = griddata(NKG2016LUabs.lat,NKG2016LUabs.lon,NKG2016LUabs.lu-NKG2016LUlev.lu,...
        TGinfo.Lat(i),  TGinfo.Lon(i), "linear")/10; %convert to cm/year

    tg{i}.dt3 =  tg{i}.dt2 + (tg{i}.y-epoch) .* TGinfo.GC(i);

end

clearvars i NAPlat epoch
%% data gap
for i=1:length(tg)
    gap(i,1)=(sum(isnan(tg{i}.dt3))/height(tg{i}))*100;
end

%% SA treatment
%  Ellipsoid correction

% t/p
a1=6378136.3;
e1=0.081819221456;
% GRS80

a2=6378137.0;
e2=0.081819190842621;

for i=1:13
    h1=(a1*(1-e1^2))./(1-(e1^2)*((sind(sa_tg{i}.lat)).^2)).^.5;
    h2=(a2*(1-e2^2))./(1-(e2^2)*((sind(sa_tg{i}.lat)).^2)).^.5;
    Ell_corr=h1-h2;

    sa_tg{i}.dt = (sa_tg{i}.dt+Ell_corr);
    clear Ell_corr
end
i
%% cut overlapped data 

for i=1:13
sa61{i}(sa61{i}.t<=max(sa6{i}.t),:) = [];
sa81{i}(sa81{i}.t<=max(sa8{i}.t),:) = [];
sa91{i}(sa91{i}.t<=max(sa9{i}.t),:) = [];
end


%% combine all sa
% load sa_cut_all

  for i=1:13
  sa1{i} = removevars(sa1{i}, "dac");
  sa1{i} = removevars(sa1{i}, "oct");
  sa2{i} = removevars(sa2{i}, "dac");
  sa2{i} = removevars(sa2{i}, "oct");
  sa3{i} = removevars(sa3{i}, "dac");
  sa3{i} = removevars(sa3{i}, "oct");
  sa4{i} = removevars(sa4{i}, "dac");
  sa4{i} = removevars(sa4{i}, "oct");
  sa5{i} = removevars(sa5{i}, "dac");
  sa5{i} = removevars(sa5{i}, "oct");
  sa6{i} = removevars(sa6{i}, "dac");
  sa6{i} = removevars(sa6{i}, "oct");
  sa7{i} = removevars(sa7{i}, "dac");
  sa7{i} = removevars(sa7{i}, "oct");
  sa8{i} = removevars(sa8{i}, "dac");
  sa8{i} = removevars(sa8{i}, "oct");
  sa9{i} = removevars(sa9{i}, "dac");
  sa9{i} = removevars(sa9{i}, "oct");
  end
%

for i=1:13
    sa_tg{i}=[sa1{i};sa2{i};sa3{i};sa4{i};sa5{i};sa6{i};sa61{i};sa7{i};sa8{i};sa9{i};sa81{i};sa91{i};sa10{i}];
    sa_tg{i} = sortrows(sa_tg{i},'t','ascend');
end

% clearvars a1 a2 e1 e2 h1 h2 sa1 sa2 sa3 sa4 sa5 sa6 sa7 sa8 sa9 sa61 sa81 sa91 sa10 Ell_corr i


%% ITRF to ETRF
load('TP.mat')

for i=1:length(tg)
    % load T/P ellipsoid
    spheroid = referenceEllipsoid('GRS 80');
    sa{i}.dt2=nan(height(sa{i}),1);

    for year=min(sa{i}.year_sa):max(sa{i}.year_sa)
        tc=year;

        % if year<2008
        %     tc=2008;
        % else
        %     tc=year;
        % end

        % convert Lat, Lon, H to X Y Z (from T/P Ellipsoid)
        [x,y,z] = geodetic2ecef(spheroid,sa{i}.lat(sa{i}.year_sa==year),sa{i}.lon(sa{i}.year_sa==year),sa{i}.dt(sa{i}.year_sa==year));


        % transfering form: ITRF2008	to: ETRF2000	 epoch: SA Year;
        X2=itrstrafo([x,y,z],'ITRF2008',tc,'ETRF2000',tc);

        % convert X Y Z to Lat Lon H (to GRS 80 Ellipsoid)
        [~,~,h] = ecef2geodetic(spheroid,X2(:,1),X2(:,2),X2(:,3));

        sa{i}.dt2(sa{i}.year_sa==year)=h*100;
        clearvars h h1 X2

    end
    clearvars year X2 x y z
end

% remove gross
for i=1:13
[m,~]=find(sa{i}.dt2<-100);
sa{i}(m,:)=[];
end
%% load final data
load('F:\Analysis\BS\trend\ASL trend\data_all_v4_13TGs_v3.mat')
% load('C:\Users\Farzad Mostafavi\OneDrive - Tallinna Tehnikaülikool\trend\ASL trend\sa_new.mat')
% load('C:\Users\Farzad Mostafavi\OneDrive - Tallinna Tehnikaülikool\trend\ASL trend\Physics.mat')

%% dt TG at SA cycle
for i=1:length(sa)
    sa{i}.date=decyear(sa{i}.t);
    for j=1:size(sa{i},1)
        if all(isnan(TG{i}.dt3(TG{i}.Time>=sa{i}.t(j)-hours(2)&TG{i}.Time<sa{i}.t(j)+hours(2))))
            sa{i}.dttg(j)=nan;
        else
            tt=decyear(TG{i}.Time(TG{i}.Time>=sa{i}.t(j)-hours(2)&TG{i}.Time<=sa{i}.t(j)+hours(2)));
            dtg=TG{i}.dt3(TG{i}.Time>=sa{i}.t(j)-hours(2)&TG{i}.Time<=sa{i}.t(j)+hours(2));
            if size(dtg(~isnan(dtg)),1)==1
                sa{i}.dttg(j)=dtg(~isnan(dtg));
            else
                sa{i}.dttg(j)=interp1(tt(~isnan(dtg)),dtg(~isnan(dtg)),sa{i}.date(j));
            end
            clearvars tt dtg
        end
    end
end

%% convert to timetable
for i=1:length(tg)
    tgs=tg{i}(tg{i}.y>1994&tg{i}.y<2023,:);
    TG{i} = retime(table2timetable(tgs), 'hourly', 'mean'); % hourly averages
    TG{i}.date = decyear(TG{i}.Time);
    TG{i}.y = year(TG{i}.Time);
    clear tgs
end

for i=1:length(sa)
    SA{i} = retime(table2timetable(sa{i}), 'hourly', 'mean'); % hourly averages
    clear tgs
end


%count data
for i=1:length(TG)
I(i,1)=size(TG{i}.obs(~isnan(TG{i}.obs)),1);
Igapped(i,1)=(size(TG{i}.obs(isnan(TG{i}.obs)),1)/size(TG{i}.obs,1))*100;


%count data 2000-2010

I(i,2)=size(TG{i}.obs(~isnan(TG{i}.obs)&TG{i}.y<2010&TG{i}.y>=2000),1);
Igapped(i,2)=(size(TG{i}.obs(isnan(TG{i}.obs)&TG{i}.y<2010&TG{i}.y>=2000),1)/size(TG{i}.obs(TG{i}.y<2010&TG{i}.y>=2000),1))*100;


%count data 2010-2020

I(i,3)=size(TG{i}.obs(~isnan(TG{i}.obs)&TG{i}.y<2020&TG{i}.y>=2010),1);
Igapped(i,3)=(size(TG{i}.obs(isnan(TG{i}.obs)&TG{i}.y<2020&TG{i}.y>=2010),1)/size(TG{i}.obs(TG{i}.y<2020&TG{i}.y>=2010),1))*100;


%count data 2000-2020

I(i,4)=size(TG{i}.obs(~isnan(TG{i}.obs)&TG{i}.y<2020&TG{i}.y>=2000),1);
Igapped(i,4)=(size(TG{i}.obs(isnan(TG{i}.obs)&TG{i}.y<2020&TG{i}.y>=2000),1)/size(TG{i}.obs(TG{i}.y<2020&TG{i}.y>=2000),1))*100;
end

for i=1:length(SA)
J(i,1)=size(SA{i}.dt(~isnan(SA{i}.dt)),1);
Jgapped(i,1)=(size(SA{i}.dt(isnan(SA{i}.dt)),1)/size(SA{i}.dt,1))*100;


%count data 2000-2010

J(i,2)=size(SA{i}.dt(~isnan(SA{i}.dt)&SA{i}.y<2010&SA{i}.y>=2000),1);
Jgapped(i,2)=(size(SA{i}.dt(isnan(SA{i}.dt)&SA{i}.y<2010&SA{i}.y>=2000),1)/size(SA{i}.dt(SA{i}.y<2010&SA{i}.y>=2000),1))*100;

%count data 2010-2020
J(i,3)=size(SA{i}.dt(~isnan(SA{i}.dt)&SA{i}.y<2020&SA{i}.y>=2010),1);
Jgapped(i,3)=(size(SA{i}.dt(isnan(SA{i}.dt)&SA{i}.y<2020&SA{i}.y>=2010),1)/size(SA{i}.dt(SA{i}.y<2020&SA{i}.y>=2010),1))*100;


%count data 2000-2020
J(i,4)=size(SA{i}.dt(~isnan(SA{i}.dt)&SA{i}.y<2020&SA{i}.y>=2000),1);
Jgapped(i,4)=(size(SA{i}.dt(isnan(SA{i}.dt)&SA{i}.y<2020&SA{i}.y>=2000),1)/size(SA{i}.dt(SA{i}.y<2020&SA{i}.y>=2000),1))*100;
end

m=unique(sa{1}.said);

for i=1:13
    % all
    JJ(i,1)=0;
    JJ(i,2)=0;
    for j=1:size(m,1)
        JJ(i,1)=JJ(i,1)+size(unique(sa{i}.id(sa{i}.said==m(j))),1);
        JJ(i,2)=JJ(i,2)+size(unique(sa{i}.id(sa{i}.said==m(j)& ~isnan(sa{i}.dttg) )),1);
    end
    %2000-2010
        JJ(i,3)=0;
        JJ(i,4)=0;
    for j=1:size(m,1)
        JJ(i,3)=JJ(i,3)+size(unique(sa{i}.id(sa{i}.said==m(j)&sa{i}.y<2010&sa{i}.y>=2000)),1);
        JJ(i,4)=JJ(i,4)+size(unique(sa{i}.id(sa{i}.said==m(j)&sa{i}.y<2010&sa{i}.y>=2000& ~isnan(sa{i}.dttg))),1);
    end
    % 2010-2020
        JJ(i,5)=0;
        JJ(i,6)=0;
    for j=1:size(m,1)
        JJ(i,5)=JJ(i,5)+size(unique(sa{i}.id(sa{i}.said==m(j)&sa{i}.y<2020&sa{i}.y>=2010)),1);
        JJ(i,6)=JJ(i,6)+size(unique(sa{i}.id(sa{i}.said==m(j)&sa{i}.y<2020&sa{i}.y>=2010& ~isnan(sa{i}.dttg))),1);
    end
    %2000-2012
    
        JJ(i,7)=0;
        JJ(i,8)=0;       
    for j=1:size(m,1)
        JJ(i,7)=JJ(i,7)+size(unique(sa{i}.id(sa{i}.said==m(j)&sa{i}.y<2020&sa{i}.y>=2000)),1);
        JJ(i,8)=JJ(i,8)+size(unique(sa{i}.id(sa{i}.said==m(j)&sa{i}.y<2020&sa{i}.y>=2000& ~isnan(sa{i}.dttg))),1);
    end
end



% count sa data
for j=1:9
for k=1:13
id=unique(sa{k}.id(sa{k}.said==j));

for i=1:length(id)
mm(i,1)=height(sa{k}.id(sa{k}.id==id(i)&sa{k}.said==1));
end

m{j}(k,1)=min(mm);
m{j}(k,2)=max(mm);
clearvars mm id

end
end

clear m
%% plot TG only
close all
for i=1:length(tg)
    figure(i)
    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=1995)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=1995),'-','color',[.7 .7 .7],'DisplayName','TG observation','LineWidth',1.5)
    hold on
    tr1=fitlm(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=1995)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=1995));
    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=1995)),tr1.Fitted,'-k','DisplayName',strcat('TG_{Linear Trend}: ',num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5)

    tr2=fitlm(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=1995&tg_mean_m{i}.ye<2000)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=1995&tg_mean_m{i}.ye<2000));
    tr3=fitlm(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2000&tg_mean_m{i}.ye<2005)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=2000&tg_mean_m{i}.ye<2005));
    tr4=fitlm(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2005&tg_mean_m{i}.ye<2010)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=2005&tg_mean_m{i}.ye<2010));
    tr5=fitlm(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2010&tg_mean_m{i}.ye<2015)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=2010&tg_mean_m{i}.ye<2015));
    tr6=fitlm(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2015&tg_mean_m{i}.ye<2020)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=2015&tg_mean_m{i}.ye<2020));


    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=1995&tg_mean_m{i}.ye<2000)),tr2.Fitted,'-b','DisplayName',strcat('TG_{95-00}: ',num2str(tr2.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5)
    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2000&tg_mean_m{i}.ye<2005)),tr3.Fitted,'-r','DisplayName',strcat('TG_{00-05}: ',num2str(tr3.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5)
    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2005&tg_mean_m{i}.ye<2010)),tr4.Fitted,'-g','DisplayName',strcat('TG_{05-10}: ',num2str(tr4.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5)
    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2010&tg_mean_m{i}.ye<2015)),tr5.Fitted,'-y','DisplayName',strcat('TG_{10-15}: ',num2str(tr5.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5)
    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2015&tg_mean_m{i}.ye<2020)),tr6.Fitted,'-m','DisplayName',strcat('TG_{15-20}: ',num2str(tr6.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5)


    xlim([min(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=1995)))-0.1,max(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye<=2022)))+0.1])

    legend show
    leg=legend;
    title(leg,TGinfo.TG(i))

    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=20; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
    set(gca,'fontname','Times New Roman','FontSize',20);

end

%% TG only rate and acceleration VLM and without VLM
% table 3

tgres=table();
% The 100(1 – α)% confidence intervals for regression coefficient
alpha=0.05;

for i=1:13
    mintime=min(sa{i}.t);
    maxtime=max(sa{i}.t);

dttg1=TG{i}.dt3(TG{i}.Time<=maxtime&TG{i}.Time>=mintime);
dttg2=TG{i}.dt3(TG{i}.Time<=maxtime&TG{i}.Time>=mintime)-((TG{i}.y(TG{i}.Time<=maxtime&TG{i}.Time>=mintime)-2000) .* TGinfo.NKG2016LU(i)');
tgti=decyear(TG{i}.Time(TG{i}.Time<=maxtime&TG{i}.Time>=mintime));



tr_tg1=fitlm(tgti,dttg1);
tr_tg2=fitlm(tgti,dttg2);

ac_tg1=fitlm(tgti.^2,dttg1); %abs
ac_tg2=fitlm(tgti.^2,dttg2); %rel

tgres.dttg1(i,1)=tr_tg1.Coefficients.Estimate(2)*10;
tgres.dttg2(i,1)=tr_tg2.Coefficients.Estimate(2)*10;

tgres.actg1(i,1)=ac_tg1.Coefficients.Estimate(2)*10;
tgres.actg2(i,1)=ac_tg2.Coefficients.Estimate(2)*10;


co_tg1=coefCI(tr_tg1,alpha);
co_tg2=coefCI(tr_tg2,alpha);

tgres.ertg1(i,1)=(co_tg1(2,2)-co_tg1(2,1))/2;
tgres.ertg2(i,1)=(co_tg2(2,2)-co_tg2(2,1))/2;

tgres.sd(i,1)=std(dttg1,'omitnan');
tgres.mean(i,1)=mean(dttg1,'omitnan');

end


%% distance to the coast

load('bs.mat')

for i=1:length(tg)
    %    find closest point to the coast of sa data
    [k,~] = dsearchn(bs,[sa{i}.lat, sa{i}.lon]);

    %measure distance to the bs coastline
    parfor j=1:length(k)
        dis(j,1)=distance([sa{i}.lat(j,1), sa{i}.lon(j,1)],[bs(k(j),1),bs(k(j),2)],referenceEllipsoid('WGS84'))/1000;
    end

    sa{i}.dist=dis;
    clearvars k dis j
end
clear bs

%% plot TG and SA trend
close all


limit=5; %select a distance from coast to remove non qualified data
col=[0 0.4470 0.7410;
    0.8500 0.3250 0.0980;   
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840;
    0 1 1;
    1 0 1
    1 .9 .8];

mission={'CS2' 'ENV' 'ER2' 'JA1' 'JA2' 'JA3' 'SRL' 'S3A' 'S3B', 'S6A'};


for i=1:length(tg)
    figure(i)


    % prepare data
    mintime=min(sa{i}.t);
    maxtime=max(sa{i}.t);
    ti=decyear(sa{i}.t(sa{i}.dist>limit));

    tg{i} = sortrows(tg{i},'Time','ascend');
    [~,ia] = unique(tg{i}.Time);
    tg{i}= tg{i}(ia,:);
    F=griddedInterpolant(decyear(tg{i}.Time),tg{i}.dt3);
    dttg=F(ti);
    dtsa=sa{i}.dt(sa{i}.dist>limit)*100;


    % DT dtrend plot
    subplot(2,3,1:3)

    %     plot instantenous DT
    plot(sa{i}.t(sa{i}.dist>limit),dtsa,'.','color',[0.3010 0.7450 0.9330],'DisplayName','SA')
    hold on
    plot(sa{i}.t(sa{i}.dist>limit),dttg,'.','color',[.7 .7 .7],'DisplayName','TG')

    %     plot linear trend

    tr2=fitlm(ti,dtsa);

    plot(sa{i}.t(sa{i}.dist>limit),tr2.Fitted,'-','color',[0 0.4470 0.7410],'DisplayName',strcat('SA_{Linear Trend}: ',...
        num2str(tr2.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5);


    tr1=fitlm(ti,dttg);
    plot(sa{i}.t(sa{i}.dist>limit),tr1.Fitted,'-k','DisplayName',strcat('TG_{Linear Trend}: ',...
        num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',1.5);


    xlim([mintime-0.5,maxtime+0.5])

    leg = legend('show','location','southeast','NumColumns',2);
    title(leg,TGinfo.TG(i))
    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; ax.FontName='Times New Roman';

    %     plot data location
    subplot(2,3,4)
    k=1;

    for j=1:10
        if ~isempty(sa{i}.lat(sa{i}.said==j))
            H(k)=geoplot(sa{i}.lat(sa{i}.said==j),sa{i}.lon(sa{i}.said==j),'.','Color',col(j,:),'DisplayName',mission{j});
            k=k+1;
        end
        hold on
    end


    geoplot(sa{i}.lat(sa{i}.dist<limit),sa{i}.lon(sa{i}.dist<limit),'.w')
    geoplot(TGinfo.Lat(i),TGinfo.Lon(i),'Marker', '^','Color','k','MarkerFaceColor','k','MarkerSize',20)


    text(TGinfo.Lat(i),TGinfo.Lon(i)+0.1, TGinfo.TG(i),'FontSize',17,'Color','k','FontWeight','Bold');
    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; ax.FontName='Times New Roman';
    legend(H(1:k-1))


    hold off


    %	multi-deg poly
    subplot(2,3,5:6)

    p1= polyfit(ti,dttg,2);
    fit1=polyval(p1,ti);
    h(1)=plot(ti,fit1,'-','color',[0 0.4470 0.7410],'LineWidth',2.5,'DisplayName',strcat('2^{nd}poly trend'));
    hold on

    p2= polyfit(ti,dttg,3);
    fit2=polyval(p2,ti);
    h(2)=plot(ti,fit2,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2.5,'DisplayName',strcat('3^{nd}poly trend'));

    p3= polyfit(ti,dttg,4);
    fit3=polyval(p3,ti);
    h(3)=plot(ti,fit3,'-','color',[0.9290 0.6940 0.1250],'LineWidth',2.5,'DisplayName',strcat('4^{th}poly trend'));

    p4= polyfit(ti,dttg,5);
    fit4=polyval(p4,ti);
    h(4)=plot(ti,fit4,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2.5,'DisplayName',strcat('5^{th}poly trend'));

    % p5= polyfit(ti,dttg,6);
    % fit5=polyval(p5,ti);
    % h(5)=plot(ti,fit5,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2.5,'DisplayName',strcat('6^{th}poly trend'));


    pp1= polyfit(ti,dtsa,2);
    fitt1=polyval(pp1,ti);
    hh(1)=plot(ti,fitt1,'--','color',[0 0.4470 0.7410],'LineWidth',2.5,'DisplayName',strcat('2^{nd}ppoly trend'));

    pp2= polyfit(ti,dtsa,3);
    fitt2=polyval(pp2,ti);
    hh(2)=plot(ti,fitt2,'--','color',[0.4660 0.6740 0.1880],'LineWidth',2.5,'DisplayName',strcat('3^{nd}ppoly trend'));

    pp3= polyfit(ti,dtsa,4);
    fitt3=polyval(pp3,ti);
    hh(3)=plot(ti,fitt3,'--','color',[0.9290 0.6940 0.1250],'LineWidth',2.5,'DisplayName',strcat('4^{thh}ppoly trend'));

    pp4= polyfit(ti,dtsa,5);
    fitt4=polyval(pp4,ti);
    hh(4)=plot(ti,fitt4,'--','color',[0.8500 0.3250 0.0980],'LineWidth',2.5,'DisplayName',strcat('5^{thh}ppoly trend'));

    % pp5= polyfit(ti,dtsa,6);
    % fitt5=polyval(pp5,ti);
    % hh(5)=plot(ti,fitt5,'--','color',[0.6350 0.0780 0.1840],'LineWidth',2.5,'DisplayName',strcat('6^{thh}ppoly trend'));

    xlim([1995-0.5,2022+0.5])
    xticks(1995:3:2023)

    %   xticks(datetime('01-Jan-1995'):calmonths(36):datetime('01-Jan-2023'))
    %   xtickformat('yy')

    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; ax.FontName='Times New Roman'; grid on

    leg=legend(h,'location','south','NumColumns',5);
    title(leg,'   — TG    --- SA')

    clearvars mintime maxtime ti dtsa dttg h p1 p2 p3 p4 p5 fit1 fit2 fit3 fit4 fit5 tr1 tr2 delta tr1 tr2 j leg H k
    clearvars pp1 pp2 pp3 pp4 pp5 fitt1 fitt2 fitt3 fitt4 fitt5 hh leg ax k ti dtsa dttg ia F

end
clearvars col limit



%% MDT monthly and yearly

% yearly MDT
% yearly mean TG
for k=1:length(TG)
    tg_mean{k}=table();
    [x1,tm] = findgroups(TG{k}.y);
    tg_mean{k}.y=tm;
    % vlm corrected
    tg_mean{k}.dttg=splitapply(@(x)mean(x,'omitnan'),TG{k}.dt3,x1); %tide system corrected
    % vlm uncorrected
    vlm=((TG{k}.y-2000).* TGinfo.NKG2016LU(k)');
    tg_mean{k}.dttgrel=splitapply(@(x)mean(x,'omitnan'),TG{k}.dt3-vlm,x1); %tide system corrected

    [m,~]=find(tg_mean{k}.y>2022);
    tg_mean{k}(m,:)=[];
    [n,~]=find(tg_mean{k}.y<1995);
    tg_mean{k}(n,:)=[];
    clearvars tm x1 m n vlm
end

% monthly mean

% TG monthly mean
%sort by date & monthly category
for k=1:length(TG)
    TG{k}.mo = month(TG{k}.Time);
    TG{k}.da = day(TG{k}.Time);
    TG{k}.ym=str2num([num2str(TG{k}.y),num2str(TG{k}.mo,'%0.2i')]);
    % tg{k}.ym2=str2num([num2str(tg{k}.y),num2str(tg{k}.mo,'%0.2i'),num2str(tg{k}.da,'%0.2i')]);
end

for k=1:length(TG)
    tg_mean_m{k}=table();
    [x1,tm] = findgroups(TG{k}.ym);
    vlm=((TG{k}.y-2000).* TGinfo.NKG2016LU(k)');

    tg_mean_m{k}.dttg=splitapply(@(x)mean(x,'omitnan'),TG{k}.dt3,x1); %tide system corrected
    tg_mean_m{k}.dttgrel=splitapply(@(x)mean(x,'omitnan'),TG{k}.dt3-vlm,x1); %relatvive

    tg_mean_m{k}.tm=tm;

    % adding back the time
    tg_mean_m{k}.ye=round((tg_mean_m{k}.tm)/100,0);
    tg_mean_m{k}.mo=round((((tg_mean_m{k}.tm)/100)-round((tg_mean_m{k}.tm)/100,0))*100,0);
    tg_mean_m{k}.time=datetime(tg_mean_m{k}.ye,tg_mean_m{k}.mo,1);


    clearvars vlm x1 tm
end
clear k


%SA monthly mean
%sort by date & monthly category
for k=1:length(SA)
    SA{k} = sortrows(SA{k},'t','ascend');
    SA{k}.y = year(SA{k}.t);
    SA{k}.mo = month(SA{k}.t);
    SA{k}.da = day(SA{k}.t);
    SA{k}.ym=str2num([num2str(SA{k}.y),num2str(SA{k}.mo,'%0.2i')]);
    % sa{k}.ym2=str2num([num2str(sa{k}.y),num2str(sa{k}.mo,'%0.2i'),num2str(sa{k}.da,'%0.2i')]);
end

for k=1:length(SA)
    sa_mean_m{k}=table();
    [x1,tm] = findgroups(SA{k}.ym);
    sa_mean_m{k}.dtsa=splitapply(@(x)mean(x,'omitnan'),SA{k}.dt2,x1); %tide system corrected
    sa_mean_m{k}.dttg=splitapply(@(x)mean(x,'omitnan'),SA{k}.dttg,x1); %tg @ SA cycle
    sa_mean_m{k}.tm=tm;
    clearvars x1 tm

    % adding back the time
    sa_mean_m{k}.ye=round((sa_mean_m{k}.tm)/100,0);
    sa_mean_m{k}.mo=round((((sa_mean_m{k}.tm)/100)-round((sa_mean_m{k}.tm)/100,0))*100,0);
    sa_mean_m{k}.time=datetime(sa_mean_m{k}.ye,sa_mean_m{k}.mo,1);
end
clear k

for k=1:13
    mean_all{k}=table();
    mean_all{k}.time=[datetime('01-Jan-1995'):calmonths(1):datetime('01-Dec-2022')]';
    for i=1:height(mean_all{k})
        if isempty(sa_mean_m{k}.dtsa(sa_mean_m{k}.time==mean_all{k}.time(i)))
            mean_all{k}.sa(i)=nan;
        else
            mean_all{k}.sa(i)=sa_mean_m{k}.dtsa(sa_mean_m{k}.time==mean_all{k}.time(i));
        end
        if isempty(tg_mean_m{k}.dttg(tg_mean_m{k}.time==mean_all{k}.time(i)))
            mean_all{k}.tg(i)=nan;
        else
            mean_all{k}.tg(i)=tg_mean_m{k}.dttg(tg_mean_m{k}.time==mean_all{k}.time(i));
        end
    end
end


%% plot gap data

% figure hourly
for i=1:length(TG)
plot(repmat(i,height(TG{i}),1),TG{i}.Time,'sb','MarkerFaceColor','b')
hold on
plot(repmat(i,(size(TG{i}.obs(isnan(TG{i}.obs)),1)),1),TG{i}.Time(isnan(TG{i}.dt3)),'s','color',[.4 .8 1],'MarkerFaceColor',[.4 .8 1])
end

for i=1:length(SA)
plot(repmat(i+0.5,height(SA{i}),1),SA{i}.t,'sr','MarkerFaceColor','r')
hold on
plot(repmat(i+0.5,(size(SA{i}.dt(isnan(SA{i}.dt)),1)),1),SA{i}.t(isnan(SA{i}.dt)),'s','color',[1 .5 .5],'MarkerFaceColor',[1 .5 .5])
end


ylim(datetime([datetime('01-Dec-1994') datetime('01-Feb-2023')]))
xticks(1:1:length(tg))
    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; ax.FontName='Times New Roman'; grid on
  
% figure monthly    
    s=10;
for i=1:length(TG)
    
plot(repmat(i,height(mean_all{i}),1),mean_all{i}.time,'sb','MarkerFaceColor','b','MarkerSize',s)
hold on
plot(repmat(i,(size(mean_all{i}.tg(isnan(mean_all{i}.tg)),1)),1),mean_all{i}.time(isnan(mean_all{i}.tg)),'s','color',[.4 .8 1],'MarkerFaceColor',[.4 .8 1],'MarkerSize',s)
end

    
for i=1:length(TG)
    
plot(repmat(i+.5,height(mean_all{i}),1),mean_all{i}.time,'s','color',[0.4660 0.6740 0.1880],	'MarkerFaceColor',[0.4660 0.6740 0.1880]	,'MarkerSize',s)
hold on
plot(repmat(i+.5,(size(mean_all{i}.sa(isnan(mean_all{i}.sa)),1)),1),mean_all{i}.time(isnan(mean_all{i}.sa)),'sr','MarkerFaceColor','r','MarkerSize',s)
end

ylim(datetime([datetime('01-Dec-1994') datetime('01-Feb-2023')]))
xticks(1:1:length(tg))
    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=20; ax.FontWeight='Bold'; ax.FontName='Times New Roman'; grid on
xtickangle(45)
xlim([0.5 14])

%% plot yearly MDT
% figure 3

close all

for k=1:length(tg)
    figure(k)
    % boxchart(TG{k}.y,TG{k}.dt2,'MarkerStyle','none','WhiskerLineColor','none')
    % VLM corrrected
    for i=1995:2022
        [m1,~]=find(TG{k}.dt3(TG{k}.y==i)==max(TG{k}.dt3(TG{k}.y==i)));
        [m2,~]=find(TG{k}.dt3(TG{k}.y==i)==min(TG{k}.dt3(TG{k}.y==i)));
        t=TG{k}.y(TG{k}.y==i);
        dt=TG{k}.dt3(TG{k}.y==i);
        plot(t,dt,'-b')
        hold on
        plot(t(m1),dt(m1),'_-b')
        plot(t(m2),dt(m2),'_-b')
    end
    h(1)=plot(tg_mean{k}.y,tg_mean{k}.dttg,'ob','MarkerFaceColor','b','MarkerSize',10,'displayname','GIA-corrected DT');

    clearvars t dt m1 m2 ax

    % VLM uncorrrected
    for i=1995:2022
        [m1,~]=find(TG{k}.dt3(TG{k}.y==i)==max(TG{k}.dt3(TG{k}.y==i)));
        [m2,~]=find(TG{k}.dt3(TG{k}.y==i)==min(TG{k}.dt3(TG{k}.y==i)));
        t=TG{k}.y(TG{k}.y==i);

        dt=TG{k}.dt3(TG{k}.y==i)-((TG{k}.y(TG{k}.y==i)-2000) .* TGinfo.NKG2016LU(k)');

        plot(t+.3,dt,'-k')
        hold on
        plot(t(m1)+.3,dt(m1),'_-k')
        plot(t(m2)+.3,dt(m2),'_-k')
    end
    h(2)=plot(tg_mean{k}.y+.3,tg_mean{k}.dttgrel,'ok','MarkerFaceColor','k','MarkerSize',10,'displayname','Relative DT');



    ylabel('Annual DT_{TG} variation [cm]','FontSize',18,'FontWeight','bold');
    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold'; grid on
    set(gca,'fontname','Times New Roman');
    % ylim([-10 10])
    xlim([1994 2023])

    xticks(1995:1:2022)
    % yline(0,'--k','LineWidth',1.5)
    xtickangle(45)
    pbaspect([1 .3 1])
    box on
    clearvars t dt m1 m2 ax
    leg=legend(h);
    title(leg,TGinfo.TG(k))
end
clear k


%% boxplot of TG absolute DT
% Figure x2
close all
for k=1:length(tg)
    figure(k)
h=boxplot(TG{k}.dt3(TG{k}.y>1994&TG{k}.y<2023),TG{k}.y(TG{k}.y>1994&TG{k}.y<2023));
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
set(gca,'fontname','Times New Roman');
set(h,{'linew'},{3})
box on
h=findobj('LineStyle','--'); set(h, 'LineStyle','-'); set(h, 'LineWidth',2); set(h, 'Color','b')
h=findobj('Marker','+'); set(h, 'Marker','none'); set(h,'MarkerFaceColor','r')
h=findobj('Marker','_'); set(h,'MarkerFaceColor','b')
pbaspect([1 .3 1])
% xticks(1995:1:2022)
% set(gca,'xtick',1995:2022);
xtickangle(45)
leg=legend(h);
title(leg,TGinfo.TG(k))
legend boxoff

end

clear k

%% decdal trend

limit=5;
res=table();
% The 100(1 – α)% confidence intervals for regression coefficient
alpha=0.05;

for i=1:length(tg)
    mintime=min(sa{i}.t);
    maxtime=max(sa{i}.t);
    
    %SA
    % decadal
    ti1=decyear(sa{i}.t(sa{i}.dist>limit&sa{i}.year_sa<2010&sa{i}.year_sa>=2000));
    ti2=decyear(sa{i}.t(sa{i}.dist>limit&sa{i}.year_sa<2020&sa{i}.year_sa>=2010));
    ti3=decyear(sa{i}.t(sa{i}.dist>limit&sa{i}.year_sa<2020&sa{i}.year_sa>=2000));
    res.n2(i,1)=size(unique(sa{i}.id(sa{i}.dist>limit&sa{i}.year_sa<2010&sa{i}.year_sa>=2000)),1);
    res.n3(i,1)=size(unique(sa{i}.id(sa{i}.dist>limit&sa{i}.year_sa<2020&sa{i}.year_sa>=2010)),1);
    
    dtsa1=sa{i}.dt2(sa{i}.dist>limit&sa{i}.year_sa<2010&sa{i}.year_sa>=2000);
    dtsa2=sa{i}.dt2(sa{i}.dist>limit&sa{i}.year_sa<2020&sa{i}.year_sa>=2010);
    dtsa3=sa{i}.dt2(sa{i}.dist>limit&sa{i}.year_sa<2020&sa{i}.year_sa>=2000);
    
    % all
    ti4=decyear(sa{i}.t(sa{i}.dist>limit));
    dtsa4=sa{i}.dt2(sa{i}.dist>limit);
    res.n1(i,1)=size(unique(sa{i}.id(sa{i}.dist>limit)),1);
    
    %TG
    dttg1=TG{i}.dt3(TG{i}.y<2010&TG{i}.y>=2000);
    dttg2=TG{i}.dt3(TG{i}.y<2020&TG{i}.y>=2010);
    dttg3=TG{i}.dt3(TG{i}.y<2020&TG{i}.y>=2000);
    dttg4=TG{i}.dt3(TG{i}.Time<=maxtime&TG{i}.Time>=mintime);
    
    %relative
    dttg11=TG{i}.dt3(TG{i}.y<2010&TG{i}.y>=2000)-((TG{i}.y(TG{i}.y<2010&TG{i}.y>=2000)-2000) .* TGinfo.NKG2016LU(i)');
    dttg22=TG{i}.dt3(TG{i}.y<2020&TG{i}.y>=2010)-((TG{i}.y(TG{i}.y<2020&TG{i}.y>=2010)-2000) .* TGinfo.NKG2016LU(i)');
    
    tgti1=TG{i}.date(TG{i}.y<2010&TG{i}.y>=2000);
    tgti2=TG{i}.date(TG{i}.y<2020&TG{i}.y>=2010);
    tgti4=TG{i}.date(TG{i}.Time<=maxtime&TG{i}.Time>=mintime);
    tgti3=TG{i}.date(TG{i}.y<2020&TG{i}.y>=2000);
    
    %TG @ SA time
    
    ti11=decyear(sa{i}.t(sa{i}.year_sa<2010&sa{i}.year_sa>=2000));
    ti12=decyear(sa{i}.t(sa{i}.year_sa<2020&sa{i}.year_sa>=2010));
    ti13=decyear(sa{i}.t(sa{i}.year_sa<2020&sa{i}.year_sa>=2000));
    
    dttg51=sa{i}.dttg(sa{i}.year_sa<2010&sa{i}.year_sa>=2000);
    dttg52=sa{i}.dttg(sa{i}.year_sa<2020&sa{i}.year_sa>=2010);
    dttg53=sa{i}.dttg(sa{i}.year_sa<2020&sa{i}.year_sa>=2000);
    dttg54=sa{i}.dttg;
    
    
    %interpolated at SA time
%     dttg51=interp1(tgti1,dttg1,ti1);
%     dttg52=interp1(tgti2,dttg2,ti2);
%     dttg53=interp1(tgti3,dttg3,ti3);
%     dttg54=interp1(tgti4,dttg4,ti4);
%     
    
    %trend
    tr_tg1=fitlm(tgti1,dttg1);
    tr_tg2=fitlm(tgti2,dttg2);
    tr_tg3=fitlm(tgti3,dttg3);
    tr_tg4=fitlm(tgti4,dttg4);
    %relative
    tr_tg11=fitlm(tgti1,dttg11);
    tr_tg22=fitlm(tgti2,dttg22);
    
    %TG @ SA time
    tr_tg51=fitlm(ti11,dttg51);
    tr_tg52=fitlm(ti12,dttg52);
    tr_tg53=fitlm(ti13,dttg53);
    tr_tg54=fitlm(sa{i}.date,dttg54);
    
    %interpolated at SA time
%     tr_tg51=fitlm(ti1,dttg51);
%     tr_tg52=fitlm(ti2,dttg52);
%     tr_tg53=fitlm(ti3,dttg53);
%     tr_tg54=fitlm(ti4,dttg54);
    
    
    tr_sa1=fitlm(ti1,dtsa1);
    tr_sa2=fitlm(ti2,dtsa2);
    tr_sa3=fitlm(ti3,dtsa3);
    tr_sa4=fitlm(ti4,dtsa4);
    
    res.dttg1(i,1)=tr_tg1.Coefficients.Estimate(2)*10;
    res.dttg2(i,1)=tr_tg2.Coefficients.Estimate(2)*10;
    res.dttg3(i,1)=tr_tg3.Coefficients.Estimate(2)*10;
    res.dttg4(i,1)=tr_tg4.Coefficients.Estimate(2)*10;
    %relative
    res.dttg11(i,1)=tr_tg11.Coefficients.Estimate(2)*10;
    res.dttg22(i,1)=tr_tg22.Coefficients.Estimate(2)*10;
    
    %TG @ SA time
    res.dttg51(i,1)=tr_tg51.Coefficients.Estimate(2)*10;
    res.dttg52(i,1)=tr_tg52.Coefficients.Estimate(2)*10;
    res.dttg53(i,1)=tr_tg53.Coefficients.Estimate(2)*10;
    res.dttg54(i,1)=tr_tg54.Coefficients.Estimate(2)*10;
    
    
    % SA
    res.dtsa1(i,1)=tr_sa1.Coefficients.Estimate(2)*10;
    res.dtsa2(i,1)=tr_sa2.Coefficients.Estimate(2)*10;
    res.dtsa3(i,1)=tr_sa3.Coefficients.Estimate(2)*10;
    res.dtsa4(i,1)=tr_sa4.Coefficients.Estimate(2)*10;
    
    
    %  coefficient standard errors
    % diag(sqrt(tr2.CoefficientCovariance))
    
    co_tg1=coefCI(tr_tg1,alpha);
    co_tg2=coefCI(tr_tg2,alpha);
    co_tg3=coefCI(tr_tg3,alpha);
    co_tg4=coefCI(tr_tg4,alpha);
    %relative
    co_tg11=coefCI(tr_tg11,alpha);
    co_tg22=coefCI(tr_tg22,alpha);
    
    co_tg51=coefCI(tr_tg51,alpha);
    co_tg52=coefCI(tr_tg52,alpha);
    co_tg53=coefCI(tr_tg53,alpha);
    co_tg54=coefCI(tr_tg54,alpha);
    
    co_sa1=coefCI(tr_sa1,alpha);
    co_sa2=coefCI(tr_sa2,alpha);
    co_sa3=coefCI(tr_sa3,alpha);
    co_sa4=coefCI(tr_sa4,alpha);
    
    res.ertg1(i,1)=(co_tg1(2,2)-co_tg1(2,1))/2;
    res.ertg2(i,1)=(co_tg2(2,2)-co_tg2(2,1))/2;
    res.ertg3(i,1)=(co_tg3(2,2)-co_tg3(2,1))/2;
    res.ertg4(i,1)=(co_tg4(2,2)-co_tg4(2,1))/2;
    
    %relative
    res.ertg11(i,1)=(co_tg11(2,2)-co_tg11(2,1))/2;
    res.ertg22(i,1)=(co_tg22(2,2)-co_tg22(2,1))/2;
    
    res.ertg51(i,1)=(co_tg51(2,2)-co_tg51(2,1))/2;
    res.ertg52(i,1)=(co_tg52(2,2)-co_tg52(2,1))/2;
    res.ertg53(i,1)=(co_tg53(2,2)-co_tg53(2,1))/2;
    res.ertg54(i,1)=(co_tg54(2,2)-co_tg54(2,1))/2;
    
    res.ersa1(i,1)=(co_sa1(2,2)-co_sa1(2,1))/2;
    res.ersa2(i,1)=(co_sa2(2,2)-co_sa2(2,1))/2;
    res.ersa3(i,1)=(co_sa3(2,2)-co_sa3(2,1))/2;
    res.ersa4(i,1)=(co_sa4(2,2)-co_sa4(2,1))/2;
    
    
    clearvars co_sa1 co_sa2 co_sa3 co_sa4 co_tg1 co_tg2 co_tg3 co_tg4 ti1 ti2 ti3  ti4 dttg1 dttg2 dttg3 dttg4 dtsa1 dtsa2 dtsa3 dtsa4
    clearvars tr_tg1 tr_tg2 tr_tg3 tr_tg4 tr_sa1 tr_sa2 tr_sa3 tr_sa4 tgti1 tgti2 tgti3 tgti4 maxtime mintime GC col  co_tg11 co_tg22 tr_tg11 tr_tg22 dttg11 dttg22
    
end

%% plot decadal trend SA vs TG
%figure 8 & Figure 9

% 1=2010-2000
% 2=2020-2010
% 3=2020-2000

% 4= all  *** for figure 9



% close all
for i=1:length(tg)
    
    figure(5)
    
    ertg=res.ertg4(i);
    dttg=res.dttg4(i);
    
    ersa=res.ersa4(i);
    dtsa=res.dtsa4(i);
    
    dttgs=res.dttg4(i);
    ertgs=res.ertg4(i);
    
%     % for figure 8
%     
%     %   % TG @SA time (light blue)
%     plot(i-0.1,[dttgs-3*(ertgs*10);dttgs+3*(ertgs*10)],'_-','color',[.4 .8 1],'LineWidth',3)
%     hold on
%     plot([i-0.1;i-0.1],[dttgs-3*(ertgs*10);dttgs+3*(ertgs*10)],'-','color',[.4 .8 1],'LineWidth',3)
%     h(1)=plot(i-0.1,dttgs,'o','color',[.4 .8 1],'MarkerFaceColor',[.4 .8 1],'MarkerSize',10,'DisplayName',strcat('TG_{@SA cycle}, RMSE=',num2str(rms(dttgs-dtsa),2),' mm/yr'));
%     
%     
%     % TG
%     plot(i,[dttg-3*(ertg*10);dttg+3*(ertg*10)],'_-b','LineWidth',3)
%     hold on
%     plot([i;i],[dttg-3*(ertg*10);dttg+3*(ertg*10)],'-b','LineWidth',3)
%     h(2)=plot(i,dttg,'ob','MarkerFaceColor','b','MarkerSize',10,'DisplayName',strcat('TG_{hourly}, RMSE=',num2str(rms(dttg-dtsa),2),' mm/yr'));
%     
%     
%     
%     % SA
%     plot(i+.1,[dtsa-3*(ersa*10);dtsa+3*((ersa*10))],'_-r','LineWidth',3)
%     hold on
%     plot([i+.1;i+.1],[dtsa-3*(ersa*10);dtsa+3*(ersa*10)],'-r','LineWidth',3)
%     h(3)=plot(i+.1,dtsa,'or','MarkerFaceColor','r','MarkerSize',10,'DisplayName','SA trend');
%     
    % for figure 9
    
     % hourly TG (dark blue)
           plot(i,[dttg-3*(ertg*10);dttg+3*(ertg*10)],'_-b','LineWidth',3)
           hold on
           plot([i;i],[dttg-3*(ertg*10);dttg+3*(ertg*10)],'-b','LineWidth',3)
           h(1)=plot(i,dttg,'ob','MarkerFaceColor','b','MarkerSize',10,'DisplayName','TG_{hourly}');
    
%        % TG @SA time (light blue)
%        plot(i,[dttgs-3*(ertgs*10);dttgs+3*(ertgs*10)],'_-','color',[.4 .8 1],'LineWidth',3)
%        hold on
%        plot([i;i],[dttgs-3*(ertgs*10);dttgs+3*(ertgs*10)],'-','color',[.4 .8 1],'LineWidth',3)
%        h(1)=plot(i,dttgs,'o','color',[.4 .8 1],'MarkerFaceColor',[.4 .8 1],'MarkerSize',10,'DisplayName','TG_{@SA cycle}');
%   
  	  %SA
      plot(i+.1,[dtsa-3*(ersa*10);dtsa+3*((ersa*10))],'_-r','LineWidth',3)
      hold on
      plot([i+.1;i+.1],[dtsa-3*(ersa*10);dtsa+3*(ersa*10)],'-r','LineWidth',3)
      h(2)=plot(i+.1,dtsa,'or','MarkerFaceColor','r','MarkerSize',10,'DisplayName','SA');
    
end

% ylabel('Trend (2000-2010) [mm/year]','FontSize',24,'FontWeight','bold');
%  ylabel('Trend (2010-2020) [mm/year]','FontSize',24,'FontWeight','bold');
%  ylabel('Trend (2000-2020) [mm/year]','FontSize',24,'FontWeight','bold');
 ylabel('Sea Level Trend [mm/year]','FontSize',24,'FontWeight','bold');

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');
grid minor
% leg=legend(h(1:3),'NumColumns', 2,'Location','SouthEast');




 leg=legend(h(1:2),'NumColumns', 2,'Location','SouthEast');
%   title(leg,strcat('RMSE=',num2str(rms(dttgs-dtsa),2),' mm/yr'))
    title(leg,strcat('RMSE=',num2str(rms(dttg-dtsa),2),' mm/yr'))


% ax.YTickLabel[];

xlim([0 length(tg)+1])
xticks(1:1:length(tg))
xticklabels(TGinfo.TG)
xtickangle(45)
xticklabels([])
% ylim([0 8])
pbaspect([1 .25 .25])

clearvars ax i


%% plot decadal trend TG @ SA time
%figure x4

%	1	 51=2010-2000
%	2	 52=2020-2010
%	3	 53=2020-2000
%	4	 54= all


close all
for i=1:length(tg)

  
ertg=res.ertg4(i);
dttg=res.dttg4(i);

ertg1=res.ertg54(i);
dttg1=res.dttg54(i);


    % hourly TG (dark blue)
     plot(i,[dttg-3*(ertg*10);dttg+3*(ertg*10)],'_-b','LineWidth',3)
     hold on
     plot([i;i],[dttg-3*(ertg*10);dttg+3*(ertg*10)],'-b','LineWidth',3)
     h(1)=plot(i,dttg,'ob','MarkerFaceColor','b','MarkerSize',10,'DisplayName','TG_{hourly}');
%     
     % TG @SA time (light blue)
     plot(i+.1,[dttg1-3*(ertg1*10);dttg1+3*(ertg1*10)],'_-','color',[.4 .8 1],'LineWidth',3)
     hold on
     plot([i+.1;i+.1],[dttg1-3*(ertg1*10);dttg1+3*(ertg1*10)],'-','color',[.4 .8 1],'LineWidth',3)
     h(2)=plot(i+.1,dttg1,'o','color',[.4 .8 1],'MarkerFaceColor',[.4 .8 1],'MarkerSize',10,'DisplayName','TG_{@SA cycle}');


end

%    ylabel('Trend (2000-2010) [mm/year]','FontSize',24,'FontWeight','bold');
%   ylabel('Trend (2010-2020) [mm/year]','FontSize',24,'FontWeight','bold');
%   ylabel('Trend (2000-2020) [mm/year]','FontSize',24,'FontWeight','bold');
 ylabel('Sea Level Trend [mm/year]','FontSize',24,'FontWeight','bold');

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');
grid minor

leg=legend(h(1:2),'NumColumns', 2,'Location','SouthEast');
title(leg,strcat('RMSE=',num2str(rms(dttg-dttg1),2),' mm/yr'))

% ax.YTickLabel[];

xlim([0 length(tg)+1])
xticks(1:1:length(tg))
xticklabels(TGinfo.TG)
xtickangle(90)
% xticklabels([])
% ylim([0 8])
pbaspect([1 .25 .25])

clearvars ax i
%% plot decadal TG trend relative vs absolute
%figure 4

% 1=2010-2000
% 2=2020-2010
% 3=2020-2000
% 4= all
% close all
figure(2)
for i=1:length(tg)
    plot(i,[res.dttg2(i)-3*(res.ertg2(i)*10);res.dttg2(i)+3*(res.ertg2(i)*10)],'_-b','LineWidth',3)
    hold on
    plot([i;i],[res.dttg2(i)-3*(res.ertg2(i)*10);res.dttg2(i)+3*(res.ertg2(i)*10)],'-b','LineWidth',3)
    plot(i,res.dttg2(i),'ob','MarkerFaceColor','b','MarkerSize',10)

    %relative
         plot(i+.1,[res.dttg22(i)-3*(res.ertg22(i)*10);res.dttg22(i)+3*(res.ertg22(i)*10)],'_-k','LineWidth',3)
         hold on
         plot([i+.1;i+.1],[res.dttg22(i)-3*(res.ertg22(i)*10);res.dttg22(i)+3*(res.ertg22(i)*10)],'-k','LineWidth',3)
         plot(i+.1,res.dttg22(i),'ok','MarkerFaceColor','k','MarkerSize',10)
% 
%     %sa
%     plot(i,[res.dtsa4(i)-3*(res.ersa4(i)*10);res.dtsa4(i)+3*(res.ersa4(i)*10)],'_-r','LineWidth',3)
%     plot([i;i],[res.dtsa4(i)-3*(res.ersa4(i)*10);res.dtsa4(i)+3*(res.ersa4(i)*10)],'-r','LineWidth',3)
%     plot(i,res.dtsa4(i),'or','MarkerFaceColor','r','MarkerSize',10)


end

 ylabel('TG trend_{2010-2020} [mm/yr]','FontSize',24,'FontWeight','bold');
% ylabel('Sea Level Trend [mm/year]','FontSize',24,'FontWeight','bold');

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');
grid minor


ylim([-7 20])
% ax.YTickLabel[];

xlim([0 length(tg)+1])
xticks(1:1:length(tg))
xticklabels(TGinfo.TG)
xtickangle(45)
% xticklabels([])

pbaspect([1 .25 .25])

clearvars ax i

%% %% plot decadal TG trend relative vs absolute
%figure 4


% 4= all
% 5= intep

close all
figure(2)
for i=1:length(tg)
    % hourly
    plot(i,[res.dttg2(i)-3*(res.ertg2(i)*10);res.dttg2(i)+3*(res.ertg2(i)*10)],'_-b','LineWidth',3)
    hold on
    plot([i;i],[res.dttg2(i)-3*(res.ertg2(i)*10);res.dttg2(i)+3*(res.ertg2(i)*10)],'-b','LineWidth',3)
    h(1)=plot(i,res.dttg2(i),'ob','MarkerFaceColor','b','MarkerSize',10,'DisplayName','TG_{hourly}');
    
    %interpolated
    plot(i+.1,[res.dttg5(i)-3*(res.ertg5(i)*10);res.dttg5(i)+3*(res.ertg5(i)*10)],'_-','color',[.4 .8 1],'LineWidth',3)
    hold on
    plot([i+.1;i+.1],[res.dttg5(i)-3*(res.ertg5(i)*10);res.dttg5(i)+3*(res.ertg5(i)*10)],'-','color',[.4 .8 1],'LineWidth',3)
    h(2)=plot(i+.1,res.dttg5(i),'o','color',[.4 .8 1],'MarkerFaceColor',[.4 .8 1],'MarkerSize',10,'DisplayName','TG_{@SA cycle}');
    
end

ylabel('TG Trend [mm/year]','FontSize',24,'FontWeight','bold');

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');
grid minor
leg=legend(h(1:2));
title(leg,strcat('RMSE=',num2str(rms(res.dttg4-res.dttg5),3),' cm'))
xlim([0 length(tg)+1])
xticks(1:1:length(tg))
xticklabels(TGinfo.TG)
xtickangle(90)
% xticklabels([])
ylim([-5 11])

pbaspect([1 .25 .25])

clearvars ax i


%% plot trend of TG (hourly vs interp)
% Figure X4


% 1=2010-2000
% 2=2020-2010
% 3=2020-2000
% 4= all
% close all
figure(2)
for i=1:length(tg)
    plot(i,[res.dttg2(i)-3*(res.ertg2(i)*10);res.dttg2(i)+3*(res.ertg2(i)*10)],'_-b','LineWidth',3)
    hold on
    plot([i;i],[res.dttg2(i)-3*(res.ertg2(i)*10);res.dttg2(i)+3*(res.ertg2(i)*10)],'-b','LineWidth',3)
    plot(i,res.dttg2(i),'ob','MarkerFaceColor','b','MarkerSize',10)

    %relative
         plot(i,[res.dttg22(i)-3*(res.ertg22(i)*10);res.dttg22(i)+3*(res.ertg22(i)*10)],'_-k','LineWidth',3)
         hold on
         plot([i;i],[res.dttg22(i)-3*(res.ertg22(i)*10);res.dttg22(i)+3*(res.ertg22(i)*10)],'-k','LineWidth',3)
         plot(i,res.dttg22(i),'ok','MarkerFaceColor','k','MarkerSize',10)
% 
%     %sa
%     plot(i,[res.dtsa4(i)-3*(res.ersa4(i)*10);res.dtsa4(i)+3*(res.ersa4(i)*10)],'_-r','LineWidth',3)
%     plot([i;i],[res.dtsa4(i)-3*(res.ersa4(i)*10);res.dtsa4(i)+3*(res.ersa4(i)*10)],'-r','LineWidth',3)
%     plot(i,res.dtsa4(i),'or','MarkerFaceColor','r','MarkerSize',10)


end

 ylabel('TG Trend (2000-2010) [mm/year]','FontSize',24,'FontWeight','bold');
% ylabel('Sea Level Trend [mm/year]','FontSize',24,'FontWeight','bold');

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');
grid minor



% ax.YTickLabel[];

xlim([0 length(tg)+1])
xticks(1:1:length(tg))
xticklabels(TGinfo.TG)
xtickangle(90)
% xticklabels([])

pbaspect([1 .25 .25])

clearvars ax i


%% categorised SA trend
% Figure 12

% mission={'CS2' 'ENV' 'ER2' 'JA1' 'JA2' 'JA3' 'SRL' 'S3A' 'S3B'};
%
% gr1=[2,7,8,9];
% gr2=[4,5,6];


limit=5;
res2=table();
% The 100(1 – α)% confidence intervals for regression coefficient
alpha=0.05;

for i=1:length(tg)

    %SA
    % categories
    ti1=decyear(SA{i}.t(SA{i}.dist>limit&SA{i}.year_sa<=2022&SA{i}.year_sa>=2003&SA{i}.said==2|SA{i}.said==7|SA{i}.said==8|SA{i}.said==9));
    ti2=decyear(SA{i}.t(SA{i}.dist>limit&SA{i}.year_sa<=2022&SA{i}.year_sa>=2003&SA{i}.said==4|SA{i}.said==5|SA{i}.said==6|SA{i}.said==5|SA{i}.said==10));

    dtsa1=SA{i}.dt2(SA{i}.dist>limit&SA{i}.year_sa<=2022&SA{i}.year_sa>=2003&SA{i}.said==2|SA{i}.said==7|SA{i}.said==8|SA{i}.said==9);
    dtsa2=SA{i}.dt2(SA{i}.dist>limit&SA{i}.year_sa<=2022&SA{i}.year_sa>=2003&SA{i}.said==4|SA{i}.said==5|SA{i}.said==6|SA{i}.said==5|SA{i}.said==10);

    %TG
    dttg1=TG{i}.dt3(TG{i}.y<=2022&TG{i}.y>=2003);
    tgti1=decyear(TG{i}.Time(TG{i}.y<=2022&TG{i}.y>=2003));


    %trend
    tr_tg1=fitlm(tgti1,dttg1);


    tr_sa1=fitlm(ti1,dtsa1);
    tr_sa2=fitlm(ti2,dtsa2);


    res2.dttg1(i,1)=tr_tg1.Coefficients.Estimate(2)*10;

    res2.dtsa1(i,1)=tr_sa1.Coefficients.Estimate(2)*10;
    res2.dtsa2(i,1)=tr_sa2.Coefficients.Estimate(2)*10;



    %  coefficient standard errors
    % diag(sqrt(tr2.CoefficientCovariance))



    co_tg1=coefCI(tr_tg1,alpha);

    co_sa1=coefCI(tr_sa1,alpha);
    co_sa2=coefCI(tr_sa2,alpha);


    res2.ertg1(i,1)=(co_tg1(2,2)-co_tg1(2,1))/2;

    res2.ersa1(i,1)=(co_sa1(2,2)-co_sa1(2,1))/2;
    res2.ersa2(i,1)=(co_sa2(2,2)-co_sa2(2,1))/2;



    clearvars co_sa1 co_sa2 co_sa3 co_tg1 co_tg2 co_tg3 ti1 ti2 ti3 dttg1 dttg2 dttg3 dtsa1 dtsa2 dtsa3 tr_tg1 tr_tg2 tr_tg3 tr_sa1 tr_sa2 tr_sa3 tgti1 tgti2 tgti3
end


%figure
for i=1:length(tg)

    %     yyaxis right

    plot(i,[res2.dttg1(i)-3*(res2.ertg1(i)*10);res2.dttg1(i)+3*(res2.ertg1(i)*10)],'_-b','LineWidth',3)
    hold on
    plot([i;i],[res2.dttg1(i)-3*(res2.ertg1(i)*10);res2.dttg1(i)+3*(res2.ertg1(i)*10)],'-b','LineWidth',3)
    plot(i,res2.dttg1(i),'ob','MarkerFaceColor','b','MarkerSize',10)

%sentinel
    plot(i+.1,[res2.dtsa1(i)-3*(res2.ersa1(i));res2.dtsa1(i)+3*((res2.ersa1(i)))],'_-r','LineWidth',3)
    hold on
    plot([i+.1;i+.1],[res2.dtsa1(i)-3*(res2.ersa1(i));res2.dtsa1(i)+3*(res2.ersa1(i))],'-r','LineWidth',3)
    plot(i+.1,res2.dtsa1(i),'or','MarkerFaceColor','r','MarkerSize',10)

%Jason
    plot(i-.1,[res2.dtsa2(i)-3*(res2.ersa2(i));res2.dtsa2(i)+3*((res2.ersa2(i)))],'_-','color',[0.9290 0.6940 0.1250],'LineWidth',3)
    hold on
    plot([i-.1;i-.1],[res2.dtsa2(i)-3*(res2.ersa2(i));res2.dtsa2(i)+3*(res2.ersa2(i))],'-','color',[0.9290 0.6940 0.1250]	,'LineWidth',3)
    plot(i-.1,res2.dtsa2(i),'o','MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerSize',10)


end

ylabel('Trend (2003-2022) [mm/year]','FontSize',24,'FontWeight','bold');

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');
% ax.YAxis(2).Color = 'k';
% ax.YTickLabel(1)=[];
% % ax.YTick(1)=[];
grid minor

ylim([0 13])
yticks(0:4:12)
xlim([0 length(tg)+1])
xticks(1:1:length(tg))
xticklabels(TGinfo.TG)
xtickangle(45)
% xticklabels([])
pbaspect([1 .25 .25])


%% results

limit=5;
res3=table();
% The 100(1 – α)% confidence intervals for regression coefficient
alpha=0.05;

for i=1:length(tg)

    % prepare data
    mintime=min(sa{i}.t);
    maxtime=max(sa{i}.t);
    ti=decyear(sa{i}.t(sa{i}.dist>limit));
    dtsa=sa{i}.dt(sa{i}.dist>limit)*100;

    F=griddedInterpolant(decyear(tg{i}.Time),tg{i}.dt3);
    dttg1=F(ti);

    dttg2=tg{i}.dt3(tg{i}.Time<=maxtime&tg{i}.Time>=mintime);
    tgti2=decyear(tg{i}.Time(tg{i}.Time<=maxtime&tg{i}.Time>=mintime));

    %trend
    tr_tg1=fitlm(ti,dttg1);
    tr_tg2=fitlm(tgti2,dttg2);
    tr_sa1=fitlm(ti,dtsa);

    res3.dttg1(i,1)=tr_tg1.Coefficients.Estimate(2)*10;
    res3.dttg2(i,1)=tr_tg2.Coefficients.Estimate(2)*10;
    res3.dtsa(i,1)=tr_sa1.Coefficients.Estimate(2)*10;

    co_tg1=coefCI(tr_tg1,alpha);
    co_tg2=coefCI(tr_tg2,alpha);
    co_sa1=coefCI(tr_sa1,alpha);

    res3.ertg1(i,1)=((co_tg1(2,2)-co_tg1(2,1))/2)*10;
    res3.ertg2(i,1)=((co_tg2(2,2)-co_tg2(2,1))/2)*10;
    res3.ersa(i,1)=((co_sa1(2,2)-co_sa1(2,1))/2)*10;

    clearvars co_tg1 co_tg2 co_tg3 tr_tg1 tr_tg2 tr_sa1 co_sa1 F mintime maxtime ti dtsa dttg1 dttg2 tgti2
end

%% envelpoing trend
% a sliding 28-year window is applied to the residuals; rate and acceleration values are estimated in
% each window using least squares to determine the range of variation of the rate and acceleration for all 28-year periods

slide=table();
limit=5;
alpha=0.05;
k=1;
for i=1:length(tg)

    for j=2000:2022
        dtsa=sa{i}.dt(sa{i}.dist>limit&sa{i}.y<=j)*100;
        tsa=decyear(sa{i}.t(sa{i}.dist>limit&sa{i}.y<=j));
        dttg=tg{i}.dt3(tg{i}.y<=j);
        ttg=decyear(tg{i}.Time(tg{i}.y<=j));

        tr1=fitlm(tsa,dtsa);
        tr2=fitlm(ttg,dttg);
        co_sa=coefCI(tr1,alpha);
        co_tg=coefCI(tr2,alpha);

        slide.sa(k,1)=tr1.Coefficients.Estimate(2)*10;
        slide.tg(k,1)=tr2.Coefficients.Estimate(2)*10;
        slide.saer(k,1)=((co_sa(2,2)-co_sa(2,1))/2)*10;
        slide.tger(k,1)=((co_tg(2,2)-co_tg(2,1))/2)*10;
        slide.tgid(k,1)=i;
        slide.y(k,1)=j;

        k=k+1;

        clearvars dtsa tsa dttg ttg tr1 tr2 co_sa co_tg
    end
end

%% sliding trend 
limit=5;
alpha=0.05;

slide2=table();
k=1;
for i=1:length(tg)
    y1=1995;
    thr2= 10; % decadal
    y2=y1+thr2;
    
    while y2<2023
        
        dtsa=sa{i}.dt2(sa{i}.dist>limit&sa{i}.y<y2&sa{i}.y>=y1);
        tsa=decyear(sa{i}.t(sa{i}.dist>limit&sa{i}.y<y2&sa{i}.y>=y1));
        dttg=tg{i}.dt3(tg{i}.y<y2&tg{i}.y>=y1);
        ttg=decyear(tg{i}.Time(tg{i}.y<y2&tg{i}.y>=y1));
        
        tr1=fitlm(tsa,dtsa);
        tr2=fitlm(ttg,dttg);
        co_sa=coefCI(tr1,alpha);
        co_tg=coefCI(tr2,alpha);
        
        slide2.sa(k,1)=tr1.Coefficients.Estimate(2)*10;
        slide2.saer(k,1)=((co_sa(2,2)-co_sa(2,1))/2)*10;
        slide2.tgid(k,1)=i;
        slide2.tgid_lat(k,1)=TGinfo.id_lat(i);
        slide2.y1(k,1)=y1;
        slide2.y2(k,1)=y2;
%         
%         % 80% data thereshold
%         if size(dttg,1)>size(dtsa,1)*0.2
            slide2.tger(k,1)=((co_tg(2,2)-co_tg(2,1))/2)*10;
            slide2.tg(k,1)=tr2.Coefficients.Estimate(2)*10;
            
%         else
%             slide2.tger(k,1)=nan;
%             slide2.tg(k,1)=nan;
%         end
%         
        
        y1=y1+1;
        y2=y1+thr2;
        k=k+1;
        
        clearvars dtsa tsa dttg ttg tr1 tr2 co_sa co_tg
        
    end
end

slide3=table();
k=1;
for i=1:length(tg)
    y1=1995;
    thr3= 5; %interdecadal
    y2=y1+thr3;
    
    while y2<2023
        
        dtsa=sa{i}.dt2(sa{i}.dist>limit&sa{i}.y<y2&sa{i}.y>=y1);
        tsa=decyear(sa{i}.t(sa{i}.dist>limit&sa{i}.y<y2&sa{i}.y>=y1));
        dttg=tg{i}.dt3(tg{i}.y<y2&tg{i}.y>=y1&tg{i}.y<y2);
        ttg=decyear(tg{i}.Time(tg{i}.y<y2&tg{i}.y>=y1));
        
        tr1=fitlm(tsa,dtsa);
        tr2=fitlm(ttg,dttg);
        co_sa=coefCI(tr1,alpha);
        co_tg=coefCI(tr2,alpha);
        
        slide3.sa(k,1)=tr1.Coefficients.Estimate(2)*10;
        slide3.saer(k,1)=((co_sa(2,2)-co_sa(2,1))/2)*10;
        slide3.tgid(k,1)=i;
        slide3.tgid_lat(k,1)=TGinfo.id_lat(i);
        slide3.y1(k,1)=y1;
        slide3.y2(k,1)=y2;
        
%          % 80% data thereshold
%         if size(dttg,1)>size(dtsa,1)*0.2
        slide3.tg(k,1)=tr2.Coefficients.Estimate(2)*10;
        slide3.tger(k,1)=((co_tg(2,2)-co_tg(2,1))/2)*10;
%         else
%            slide3.tg(k,1)=nan;
%         slide3.tger(k,1)=nan;
%         end
%         
        y1=y1+1;
        y2=y1+thr3;
        k=k+1;
        
        clearvars dtsa tsa dttg ttg tr1 tr2 co_sa co_tg
        
    end
end


slide4=table();
k=1;
for i=1:length(tg)
    y1=1995;
    thr3= 7; %interdecadal
    y2=y1+thr3;
    
    while y2<2023
        
        dtsa=sa{i}.dt2(sa{i}.dist>limit&sa{i}.y<y2&sa{i}.y>=y1);
        tsa=decyear(sa{i}.t(sa{i}.dist>limit&sa{i}.y<y2&sa{i}.y>=y1));
        dttg=tg{i}.dt3(tg{i}.y<y2&tg{i}.y>=y1&tg{i}.y<y2);
        ttg=decyear(tg{i}.Time(tg{i}.y<y2&tg{i}.y>=y1));
        
        tr1=fitlm(tsa,dtsa);
        tr2=fitlm(ttg,dttg);
        co_sa=coefCI(tr1,alpha);
        co_tg=coefCI(tr2,alpha);
        
        slide4.sa(k,1)=tr1.Coefficients.Estimate(2)*10;
        slide4.saer(k,1)=((co_sa(2,2)-co_sa(2,1))/2)*10;
        slide4.tgid(k,1)=i;
        slide4.tgid_lat(k,1)=TGinfo.id_lat(i);
        slide4.y1(k,1)=y1;
        slide4.y2(k,1)=y2;
        
%          % 80% data thereshold
%         if size(dttg,1)>size(dtsa,1)*0.2
        slide4.tg(k,1)=tr2.Coefficients.Estimate(2)*10;
        slide4.tger(k,1)=((co_tg(2,2)-co_tg(2,1))/2)*10;
%         else
%       slide4.tg(k,1)=nan;
%       slide4.tger(k,1)=nan;
%         end
%         
        y1=y1+1;
        y2=y1+thr3;
        k=k+1;
        
        clearvars dtsa tsa dttg ttg tr1 tr2 co_sa co_tg
        
    end
end

slide4=table();
k=1;
for i=1:length(tg)
    y1=1995;
    thr3= 7; %interdecadal
    y2=y1+thr3;
    
    while y2<2023
        
        dtsa=sa{i}.dt2(sa{i}.dist>limit&sa{i}.y<y2&sa{i}.y>=y1);
        tsa=decyear(sa{i}.t(sa{i}.dist>limit&sa{i}.y<y2&sa{i}.y>=y1));
        dttg=tg{i}.dt3(tg{i}.y<y2&tg{i}.y>=y1&tg{i}.y<y2);
        ttg=decyear(tg{i}.Time(tg{i}.y<y2&tg{i}.y>=y1));
        
        tr1=fitlm(tsa,dtsa);
        tr2=fitlm(ttg,dttg);
        co_sa=coefCI(tr1,alpha);
        co_tg=coefCI(tr2,alpha);
        
        slide4.sa(k,1)=tr1.Coefficients.Estimate(2)*10;
        slide4.saer(k,1)=((co_sa(2,2)-co_sa(2,1))/2)*10;
        slide4.tgid(k,1)=i;
        slide4.tgid_lat(k,1)=TGinfo.id_lat(i);
        slide4.y1(k,1)=y1;
        slide4.y2(k,1)=y2;
        
%          % 80% data thereshold
%         if size(dttg,1)>size(dtsa,1)*0.2
        slide4.tg(k,1)=tr2.Coefficients.Estimate(2)*10;
        slide4.tger(k,1)=((co_tg(2,2)-co_tg(2,1))/2)*10;
%         else
%       slide4.tg(k,1)=nan;
%       slide4.tger(k,1)=nan;
%         end
%         
        y1=y1+1;
        y2=y1+thr3;
        k=k+1;
        
        clearvars dtsa tsa dttg ttg tr1 tr2 co_sa co_tg
        
    end
end

clearvars k i j alpha limit
clc


% y=10

j=1;
for i=1995:2012
slide2_mean(j,1)=mean(slide2.sa(slide2.y1==i));
slide2_mean(j,2)=mean(slide2.tg(slide2.y1==i));
j=j+1;
end
% y=5

j=1;
for i=1995:2017
slide3_mean(j,1)=mean(slide3.sa(slide3.y1==i));
slide3_mean(j,2)=mean(slide3.tg(slide3.y1==i));
j=j+1;
end

% y=7
j=1;
for i=1995:2015
slide4_mean(j,1)=mean(slide4.sa(slide4.y1==i));
slide4_mean(j,2)=mean(slide4.tg(slide4.y1==i));
j=j+1;
end

clear j
%% sliding physics 

slide_ph1=table();
k=1;
    y1=1995;
    thr2= 10; % decadal
    y2=y1+thr2;
    
    while y2<2022
        
        so=ph.so(ph.y<y2&ph.y>=y1);
        th=ph.thetao(ph.y<y2&ph.y>=y1);
        nao=ph.nao(ph.y<y2&ph.y>=y1);
        uo=ph.uo(ph.y<y2&ph.y>=y1);
        vo=ph.vo(ph.y<y2&ph.y>=y1);

        t=ph.y(ph.y<y2&ph.y>=y1);
       
        
        tr1=fitlm(t,so);
        tr2=fitlm(t,th);
        tr3=fitlm(t,nao);
        tr4=fitlm(t,uo);
        tr5=fitlm(t,vo);
        
        slide_ph1.so(k,1)=tr1.Coefficients.Estimate(2)*10;
        slide_ph1.th(k,1)=tr2.Coefficients.Estimate(2)*10;
        slide_ph1.nao(k,1)=tr3.Coefficients.Estimate(2)*10;
        slide_ph1.uo(k,1)=tr4.Coefficients.Estimate(2)*10;
        slide_ph1.vo(k,1)=tr5.Coefficients.Estimate(2)*10;

        slide_ph1.y1(k,1)=y1;
        slide_ph1.y2(k,1)=y2;
     
        y1=y1+1;
        y2=y1+thr2;
        k=k+1;
        
        clearvars so th nao t tr1 tr2 tr3
        
    end


slide_ph2=table();
k=1;
    y1=1995;
    thr3= 5; %interdecadal
    y2=y1+thr3;
    
    while y2<2022
        
        so=ph.so(ph.y<y2&ph.y>=y1);
        th=ph.thetao(ph.y<y2&ph.y>=y1);
        nao=ph.nao(ph.y<y2&ph.y>=y1);
        uo=ph.uo(ph.y<y2&ph.y>=y1);
        vo=ph.vo(ph.y<y2&ph.y>=y1);

        t=ph.y(ph.y<y2&ph.y>=y1);
       
        
        tr1=fitlm(t,so);
        tr2=fitlm(t,th);
        tr3=fitlm(t,nao);
        tr4=fitlm(t,uo);
        tr5=fitlm(t,vo);
        
        slide_ph2.so(k,1)=tr1.Coefficients.Estimate(2)*10;
        slide_ph2.th(k,1)=tr2.Coefficients.Estimate(2)*10;
        slide_ph2.nao(k,1)=tr3.Coefficients.Estimate(2)*10;
        slide_ph2.uo(k,1)=tr4.Coefficients.Estimate(2)*10;
        slide_ph2.vo(k,1)=tr5.Coefficients.Estimate(2)*10;


        slide_ph2.y1(k,1)=y1;
        slide_ph2.y2(k,1)=y2;
        
        
        y1=y1+1;
        y2=y1+thr3;
        k=k+1;
        
        clearvars so th nao t tr1 tr2 tr3
        
    end


slide_ph3=table();
k=1;
    y1=1995;
    thr3= 7; %interdecadal
    y2=y1+thr3;
    
    while y2<2022
        
        so=ph.so(ph.y<y2&ph.y>=y1);
        th=ph.thetao(ph.y<y2&ph.y>=y1);
        nao=ph.nao(ph.y<y2&ph.y>=y1);
        uo=ph.uo(ph.y<y2&ph.y>=y1);
        vo=ph.vo(ph.y<y2&ph.y>=y1);

        t=ph.y(ph.y<y2&ph.y>=y1);
       
        
        tr1=fitlm(t,so);
        tr2=fitlm(t,th);
        tr3=fitlm(t,nao);
        tr4=fitlm(t,uo);
        tr5=fitlm(t,vo);
        
        slide_ph3.so(k,1)=tr1.Coefficients.Estimate(2)*10;
        slide_ph3.th(k,1)=tr2.Coefficients.Estimate(2)*10;
        slide_ph3.nao(k,1)=tr3.Coefficients.Estimate(2)*10;
        slide_ph3.uo(k,1)=tr4.Coefficients.Estimate(2)*10;
        slide_ph3.vo(k,1)=tr5.Coefficients.Estimate(2)*10;


        slide_ph3.y1(k,1)=y1;
        slide_ph3.y2(k,1)=y2;
        
        
        y1=y1+1;
        y2=y1+thr3;
        k=k+1;
        
        clearvars so th nao t tr1 tr2 tr3
        
    end


clearvars k i j alpha limit thr2 thr3
clc

%% plot enveloping trend (from 2014)
%figure 10
close all

for i=1:length(tg)
    figure(i)
    plot(slide.y(slide.tgid==i),slide.sa(slide.tgid==i),'-r','linewidth',1.5)
    hold on
    plot(slide.y(slide.tgid==i),slide.sa(slide.tgid==i)+slide.saer(slide.tgid==i),'--r')
    plot(slide.y(slide.tgid==i),slide.sa(slide.tgid==i)-slide.saer(slide.tgid==i),'--r')


    plot(slide.y(slide.tgid==i),slide.tg(slide.tgid==i),'-b','linewidth',1.5)
    hold on
    plot(slide.y(slide.tgid==i),slide.tg(slide.tgid==i)+slide.tger(slide.tgid==i),'--b')
    plot(slide.y(slide.tgid==i),slide.tg(slide.tgid==i)-slide.tger(slide.tgid==i),'--b')


    xlim([2014 2023.5])
    ylabel('Trend [mm/year]','FontSize',24,'FontWeight','bold');

    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; ax.FontName='Times New Roman';
    xticks(2010:1:2022)
    if i~=12
    xticklabels([])
    end
    
    box off
    grid on
    pbaspect([1 .25 .25])

end


%% plot enveloping 2

for i=7:length(tg)
    
    plot(slide.y(slide.tgid==i),slide.tg(slide.tgid==i),'-','linewidth',1.5,'DisplayName',strcat('TG_{id}',num2str(i)))
    hold on
end

xlim([2000 2022])
ylabel('Enveloping Trend_{TG} [mm/year]','FontSize',24,'FontWeight','bold');

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; ax.FontName='Times New Roman';
xticks(2000:2:2022)

    box on
    grid minor
    grid on
    pbaspect([1 .3 .5])



%% plot sliding

% decadal

for i=1:length(tg)
scatter(slide2.tgid(slide2.tgid==i),slide2.y1(slide2.tgid==i),800,slide2.tg(slide2.tgid==i),"square",'filled')
hold on
end
c=colorbar;
%   c.Label.String = 'SA sliding trend_{decadal} [mm/year]';
c.Label.String = 'TG sliding trend_{decadal} [mm/year]';

colormap(crameri('vik'))
xticks(1:1:length(tg))
yticks(min(slide2.y1(slide2.tgid==i)):1:max(slide2.y1(slide2.tgid==i)));
yticklabels(deca)
caxis([-15 15])
xlim([0 14])
ylim([1994 2013])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.5 1 1])




%% plot sliding
close all

figure(2)
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
ylim([1994 2013])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])



%% plot sliding


% inter-decadal

for i=1:length(tg)
scatter(slide3.tgid(slide3.tgid==i),slide3.y1(slide3.tgid==i),800,slide3.tg(slide3.tgid==i),"square",'filled')
hold on
end
c=colorbar;
% c.Label.String = 'SA sliding trend_{inter-decadal} [mm/year]';
c.Label.String = 'TG sliding trend_{inter-decadal} [mm/year]';

colormap(crameri('vik'))
xticks(1:1:length(tg))
yticks(min(slide3.y1(slide3.tgid==i)):1:max(slide3.y1(slide3.tgid==i)));
yticklabels(interd)
caxis([-45 70])
xlim([0 14])
ylim([1994 2018])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.5 1 1])



figure(2)
scatter(repmat(1,1,23),slide3.y1(slide3.tgid==1),800,slide3_mean(:,1),"square",'filled')
hold on
scatter(repmat(2,1,23),slide3.y1(slide3.tgid==1),800,slide3_mean(:,2),"square",'filled')
c=colorbar;
c.Label.String = 'mean sliding trend_{decadal} [mm/year]';

colormap(crameri('vik'))
xticks(1:2)
xticklabels(['SA';'TG'])
yticks(min(slide3.y1(slide3.tgid==1)):1:max(slide3.y1(slide3.tgid==1)));
yticklabels(interd)
caxis([-45 70])
xlim([0 3])
ylim([1994 2018])
box on
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
pbaspect([.2 1 1])


%% plot enveloping trend error @ TG locations map
%figure 10

limit=5;
latlim = [53 67];
lonlim = [10 31];
ax = usamap(latlim, lonlim);
hold on
geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
setm(ax, 'FFaceColor', [1 1 1])
setm(gca,'FLineWidth',5,'Grid','on','FontSize',18,'fontweight','bold')

for i=1:length(tg)

    rmse=rms(slide.sa(slide.tgid==i)-slide.tg(slide.tgid==i));
    r=corrcoef(slide.sa(slide.tgid==i),slide.tg(slide.tgid==i));
    rsquare=r(1,2)*100;

    scatterm(TGinfo.Lat(i), TGinfo.Lon(i),rmse*200,rsquare,'filled','o','MarkerEdgeColor','k','LineWidth',1.5);
    hold on
    clearvars rmse rsquare r
end

scatterm(55.8,29,5*200,.9,'o','MarkerEdgeColor','k','LineWidth',2);

colormap(crameri('hawaii'))

c=colorbar;
c.Label.String = 'R^2 (SA, TG trend)';


set(c,'position',[.71 .2 .015 .7],'FontSize',18,'FontWeight','bold')

% colormap(flipud(hot))

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';
xLoc =1.3001e+05;
yLoc =7.0876e+06;
scaleruler('Units', 'km', 'RulerStyle', 'patches', ...
    'XLoc', xLoc, 'YLoc', yLoc,"FontSize",10);

clearvars n1 x y ax tg_selected

%% compare gridded products and instantuneous dt

% compare with res 4th coloumn

%load gridded products
% load('F:\Analysis\BS\trend\data\gridded data\gridded.mat')

limit=5;
alpha=0.05;
res4=table();
t = datetime('2019-06-01','InputFormat','yyyy-MM-dd');

for i=1:length(tg)
    %gridded global altimetry
    tr_cp=fitlm(decyear(tg_cp{i}.t(tg_cp{i}.t<=t)),(tg_cp{i}.sla(tg_cp{i}.t<=t))*100);
    res4.trcp(i,1)=tr_cp.Coefficients.Estimate(2)*10;
    co_cp=coefCI(tr_cp,alpha);
    res4.ercp(i,1)=((co_cp(2,2)-co_cp(2,1))/2)*10;
    %gridded baltic nemo
    tr_cp2=fitlm(decyear(tg_cp2{i}.t(tg_cp2{i}.t<=t)),(tg_cp2{i}.sla(tg_cp2{i}.t<=t))*100);
    res4.trcp2(i,1)=tr_cp2.Coefficients.Estimate(2)*10;
    co_cp2=coefCI(tr_cp2,alpha);
    res4.ercp2(i,1)=((co_cp2(2,2)-co_cp2(2,1))/2)*10;
    %gridded baltic seal
    tr_se=fitlm(decyear(tg_se{i}.t),((tg_se{i}.ssh-tg_se{i}.N))*100);
    res4.trse(i,1)=tr_se.Coefficients.Estimate(2)*10;
    co_se=coefCI(tr_se,alpha);
    res4.erse(i,1)=((co_se(2,2)-co_se(2,1))/2)*10;
    % tg obs
    tr_tg=fitlm(decyear(tg{i}.Time(tg{i}.Time<=t)),tg{i}.dt3(tg{i}.Time<=t));
    res4.trtg(i,1)=tr_tg.Coefficients.Estimate(2)*10;
    co_tg=coefCI(tr_tg,alpha);
    res4.ertg(i,1)=((co_tg(2,2)-co_tg(2,1))/2)*10;
    % sa mmeasurment
    tr_sa=fitlm(decyear(sa{i}.t(sa{i}.dist>limit&sa{i}.t<=t)),(sa{i}.dt(sa{i}.dist>limit&sa{i}.t<=t))*100);
    res4.trsa(i,1)=tr_sa.Coefficients.Estimate(2)*10;
    co_sa=coefCI(tr_sa,alpha);
    res4.ersa(i,1)=((co_sa(2,2)-co_sa(2,1))/2)*10;

    clearvars tr_se co_se tr_cp co_cp tr_cp2 co_cp2
end

clearvars t


%% plot comparision (gridded vs instantenous)
%figure 13

rmsese=rms(res4.trse-res4.trtg,'omitnan');
rmsecp1=rms(res4.trcp-res4.trtg,'omitnan');
rmsecp2=rms(res4.trcp2-res4.trtg,'omitnan');
rmsesa=rms(res4.trsa-res4.trtg,'omitnan');

% remove nan values from res4 for comparision
% res5=res4;
% res5([2:3,6,8:10],:) = [];
% rmsese=rms(res5.trse-res5.trtg,'omitnan');
% rmsecp1=rms(res5.trcp-res5.trtg,'omitnan');
% rmsecp2=rms(res5.trcp2-res5.trtg,'omitnan');
% rmsesa=rms(res5.trsa-res5.trtg,'omitnan');


for i=1:length(tg)
    %tg
    plot(i,[res4.trtg(i)-3*(res4.ertg(i));res4.trtg(i)+3*(res4.ertg(i))],'_-b','LineWidth',3)
    hold on
    if i==length(tg)
        h(5)=plot([i;i],[res4.trtg(i)-3*(res4.ertg(i));res4.trtg(i)+3*(res4.ertg(i))],'-b','LineWidth',3, 'DisplayName','TG trend');
    else
        plot([i;i],[res4.trtg(i)-3*(res4.ertg(i));res4.trtg(i)+3*(res4.ertg(i))],'-b','LineWidth',3)
    end
    plot(i,res4.trtg(i),'ob','MarkerFaceColor','b','MarkerSize',10)

    %sa
    plot(i+.1,[res4.trsa(i)-3*(res4.ersa(i));res4.trsa(i)+3*(res4.ersa(i))],'_-r','LineWidth',3)
    if i==length(tg)
        h(1)=plot([i+.1;i+.1],[res4.trsa(i)-3*(res4.ersa(i));res4.trsa(i)+3*(res4.ersa(i))],'-r','LineWidth',3, 'DisplayName', strcat('SA trend. RMSE=',num2str(rmsese),'[mm/y]'));
    else
        plot([i+.1;i+.1],[res4.trsa(i)-3*(res4.ersa(i));res4.trsa(i)+3*(res4.ersa(i))],'-r','LineWidth',3)
    end
    plot(i+.1,res4.trsa(i),'or','MarkerFaceColor','r','MarkerSize',10)

    %cp1
    plot(i-.1,[res4.trcp(i)-3*(res4.ercp(i));res4.trcp(i)+3*(res4.ercp(i))],'_-','color',[0.4940 0.1840 0.5560],'LineWidth',3)
    if i==length(tg)
        h(2)=plot([i-.1;i-.1],[res4.trcp(i)-3*(res4.ercp(i));res4.trcp(i)+3*(res4.ercp(i))],'-','color',[0.4940 0.1840 0.5560],'LineWidth',3, 'DisplayName', strcat('CP1 trend. RMSE=',num2str(rmsecp1),'[mm/y]'));
    else
        plot([i-.1;i-.1],[res4.trcp(i)-3*(res4.ercp(i));res4.trcp(i)+3*(res4.ercp(i))],'-','color',[0.4940 0.1840 0.5560],'LineWidth',3)
    end
    plot(i-.1,res4.trcp(i),'o','color',[0.4940 0.1840 0.5560],'MarkerFaceColor',[0.4940 0.1840 0.5560],'MarkerSize',10)

    % cp2
    plot(i-.2,[res4.trcp2(i)-3*(res4.ercp2(i));res4.trcp2(i)+3*(res4.ercp2(i))],'_-','color',[0.4660 0.6740 0.1880],'LineWidth',3)
    if i==length(tg)
        h(3)=plot([i-.2;i-.2],[res4.trcp2(i)-3*(res4.ercp2(i));res4.trcp2(i)+3*(res4.ercp2(i))],'-','color',[0.4660 0.6740 0.1880],'LineWidth',3, 'DisplayName', strcat('CP2 trend. RMSE=',num2str(rmsecp2),'[mm/y]'));
    else
        plot([i-.2;i-.2],[res4.trcp2(i)-3*(res4.ercp2(i));res4.trcp2(i)+3*(res4.ercp2(i))],'-','color',[0.4660 0.6740 0.1880],'LineWidth',3)
    end
    plot(i-.2,res4.trcp2(i),'o','color',[0.4660 0.6740 0.1880],'MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerSize',10)

    % se
    plot(i+.2,[res4.trse(i)-3*(res4.erse(i));res4.trse(i)+3*(res4.erse(i))],'_-','color',[0.9290 0.6940 0.1250],'LineWidth',3)
    if i==length(tg)
        h(4)=plot([i+.2;i+.2],[res4.trse(i)-3*(res4.erse(i));res4.trse(i)+3*(res4.erse(i))],'-','color',[0.9290 0.6940 0.1250],'LineWidth',3, 'DisplayName', strcat('SE trend. RMSE=',num2str(rmsese),'[mm/y]'));
    else
        plot([i+.2;i+.2],[res4.trse(i)-3*(res4.erse(i));res4.trse(i)+3*(res4.erse(i))],'-','color',[0.9290 0.6940 0.1250],'LineWidth',3)
    end
    plot(i+.2,res4.trse(i),'o','color',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerSize',10)

end

for i=1:length(tg)

    ertg=res.ertg4(i);
    dttg=res.dttg4(i);

    ersa=res.ersa4(i);
    dtsa=res.dtsa4(i);


    % hourly TG (dark blue)
    % plot(i+.3,[dttg-3*(ertg*10);dttg+3*(ertg*10)],'_-','LineWidth',3,'color',[0 0 1 .2])
    % hold on
    % plot([i+.3;i+.3],[dttg-3*(ertg*10);dttg+3*(ertg*10)],'-','LineWidth',3,'color',[0 0 1 .2])
    h(1)=plot(i+.3,dttg,'^b','MarkerFaceColor','b','MarkerSize',10,'DisplayName','TG_{1995-2022}');


    %SA
    % plot(i-.3,[dtsa-3*(ersa*10);dtsa+3*((ersa*10))],'_-','LineWidth',3,'color',[1 0 0 .2])
    % hold on
    % plot([i-.3;i-.3],[dtsa-3*(ersa*10);dtsa+3*(ersa*10)],'-','LineWidth',3,'color',[1 0 0 .2])
    h(2)=plot(i-.3,dtsa,'^r','MarkerFaceColor','r','MarkerSize',10,'DisplayName','SA_{1995-2022}');

end


ylabel('Trend [mm/year]','FontSize',24,'FontWeight','bold');

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; grid on
set(gca,'fontname','Times New Roman');
% ax.YAxis(2).Color = 'k';
% ax.YTickLabel(1)=[];
% % ax.YTick(1)=[];
% grid minor

legend(h(1:5))
xlim([0 length(tg)+1])
xticks(1:1:length(tg))
xticklabels(TGinfo.TG)
xtickangle(45)
% xticklabels([])
pbaspect([1 .25 .25])


clearvars h ax rmsese rmsecp1 rmsecp2 rmsesa

%% TG trend and accelaration (abs vs rel)
alpha=0.05;
tgres=table();
for i=1:length(tg)

    mintime=min(sa{i}.t);
    maxtime=max(sa{i}.t);

    x=decyear(tg{i}.Time(tg{i}.Time<=maxtime&tg{i}.Time>=mintime));
    y1=tg{i}.dt3(tg{i}.Time<=maxtime&tg{i}.Time>=mintime);
    y2=tg{i}.dt3(tg{i}.Time<=maxtime&tg{i}.Time>=mintime)-((tg{i}.y(tg{i}.Time<=maxtime&tg{i}.Time>=mintime)-2000) .* TGinfo.NKG2016LU(i)');
    %trend
    tr1=fitlm(x,y1);
    tr2=fitlm(x,y2);

    tgres.dttg1(i,1)=tr1.Coefficients.Estimate(2)*10;
    co_tg1=coefCI(tr1,alpha);
    tgres.ertg1(i,1)=((co_tg1(2,2)-co_tg1(2,1))/2)*10;

    tgres.dttg2(i,1)=tr2.Coefficients.Estimate(2)*10;
    co_tg2=coefCI(tr2,alpha);
    tgres.ertg2(i,1)=((co_tg2(2,2)-co_tg2(2,1))/2)*10;


    %accelaration
    tr3=fitlm(x.^2,y1);
    tr4=fitlm(x.^2,y2);

    tgres.actg1(i,1)=tr3.Coefficients.Estimate(2)*10;
    co_actg1=coefCI(tr3,alpha);
    tgres.acertg1(i,1)=((co_actg1(2,2)-co_actg1(2,1))/2)*10;

    tgres.actg2(i,1)=tr4.Coefficients.Estimate(2)*10;
    co_actg2=coefCI(tr4,alpha);
    tgres.acertg2(i,1)=((co_actg2(2,2)-co_actg2(2,1))/2)*10;

    clearvars x y1 y2 tr1 tr2 co_tg1 co_tg2
end
%% annual and inter-decadal trend (13 mont moving avarage)
% figure 6
close all
marates=table();
for i=1:length(tg)

    figure(i)
    % 13 month moving avarage
    T1 = length(sa_mean_m{i}.dtsa);
    sW13 = [1/24; repmat(1/12,11,1); 1/24];
    yS = conv(sa_mean_m{i}.dtsa,sW13,'same');
    yS(1:6) = yS(7); yS(T1-5:T1) = yS(T1-6);

    T2 = length(tg_mean_m{i}.dttg);
    yT = conv(tg_mean_m{i}.dttg,sW13,'same');
    yT(1:6) = yT(7); yT(T2-5:T2) = yT(T2-6);

    trtg=fitlm(decyear(tg_mean_m{i}.time),tg_mean_m{i}.dttg);
    trsa=fitlm(decyear(sa_mean_m{i}.time),sa_mean_m{i}.dtsa);


    h(1)=plot(decyear(sa_mean_m{i}.time),yS,'-r','LineWidth',1.5,'DisplayName','SA');
    hold on
    h(2)=plot(decyear(tg_mean_m{i}.time),yT,'-b','LineWidth',1.5,'DisplayName','TG');


    % h(1)=plot(decyear(sa_mean_m{i}.time),yS,'-r','LineWidth',1.5,'DisplayName',strcat('SA trend= ',num2str(trsa.Coefficients.Estimate(2)*10,2),' mm/year'));
    % hold on
    % h(2)=plot(decyear(tg_mean_m{i}.time),yT,'-b','LineWidth',1.5,'DisplayName',strcat('TG trend=',num2str(trtg.Coefficients.Estimate(2)*10,2),' mm/year'));

    % plot(decyear(tg_mean_m{i}.time),trtg.Fitted,'-','color',[0 0.4470 0.7410],'LineWidth',1)
    % plot(decyear(sa_mean_m{i}.time),trsa.Fitted,'-','color',[0.6350 0.0780 0.1840],'LineWidth',1)



    tr1=fitlm(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=1995&tg_mean_m{i}.ye<2000)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=1995&tg_mean_m{i}.ye<2000));
    tr2=fitlm(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2000&tg_mean_m{i}.ye<2005)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=2000&tg_mean_m{i}.ye<2005));
    tr3=fitlm(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2005&tg_mean_m{i}.ye<2010)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=2005&tg_mean_m{i}.ye<2010));
    tr4=fitlm(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2010&tg_mean_m{i}.ye<2015)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=2010&tg_mean_m{i}.ye<2015));
    tr5=fitlm(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2015&tg_mean_m{i}.ye<2020)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=2015&tg_mean_m{i}.ye<2020));

    tr11=fitlm(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=1995&sa_mean_m{i}.ye<2000)),sa_mean_m{i}.dtsa(sa_mean_m{i}.ye>=1995&sa_mean_m{i}.ye<2000));
    tr12=fitlm(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=2000&sa_mean_m{i}.ye<2005)),sa_mean_m{i}.dtsa(sa_mean_m{i}.ye>=2000&sa_mean_m{i}.ye<2005));
    tr13=fitlm(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=2005&sa_mean_m{i}.ye<2010)),sa_mean_m{i}.dtsa(sa_mean_m{i}.ye>=2005&sa_mean_m{i}.ye<2010));
    tr14=fitlm(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=2010&sa_mean_m{i}.ye<2015)),sa_mean_m{i}.dtsa(sa_mean_m{i}.ye>=2010&sa_mean_m{i}.ye<2015));
    tr15=fitlm(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=2015&sa_mean_m{i}.ye<2020)),sa_mean_m{i}.dtsa(sa_mean_m{i}.ye>=2015&sa_mean_m{i}.ye<2020));


    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=1995&tg_mean_m{i}.ye<2000)),tr1.Fitted,'-','color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2000&tg_mean_m{i}.ye<2005)),tr2.Fitted,'-','color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2005&tg_mean_m{i}.ye<2010)),tr3.Fitted,'-','color',[0.4940 0.1840 0.5560],'LineWidth',1.5)
    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2010&tg_mean_m{i}.ye<2015)),tr4.Fitted,'-','color',[0.4660 0.6740 0.1880],'LineWidth',1.5)
    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2015&tg_mean_m{i}.ye<2020)),tr5.Fitted,'-','color',[0.3010 0.7450 0.9330],'LineWidth',1.5)


    plot(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=1995&sa_mean_m{i}.ye<2000)),tr11.Fitted,'--','color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
    plot(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=2000&sa_mean_m{i}.ye<2005)),tr12.Fitted,'--','color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
    plot(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=2005&sa_mean_m{i}.ye<2010)),tr13.Fitted,'--','color',[0.4940 0.1840 0.5560],'LineWidth',1.5)
    plot(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=2010&sa_mean_m{i}.ye<2015)),tr14.Fitted,'--','color',[0.4660 0.6740 0.1880],'LineWidth',1.5)
    plot(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=2015&sa_mean_m{i}.ye<2020)),tr15.Fitted,'--','color',[0.3010 0.7450 0.9330],'LineWidth',1.5)

    % TG
    marates.tr1(i)=tr1.Coefficients.Estimate(2)*10;
    marates.tr2(i)=tr2.Coefficients.Estimate(2)*10;
    marates.tr3(i)=tr3.Coefficients.Estimate(2)*10;
    marates.tr4(i)=tr4.Coefficients.Estimate(2)*10;
    marates.tr5(i)=tr5.Coefficients.Estimate(2)*10;
    % SA
    marates.tr11(i)=tr11.Coefficients.Estimate(2)*10;
    marates.tr12(i)=tr12.Coefficients.Estimate(2)*10;
    marates.tr13(i)=tr13.Coefficients.Estimate(2)*10;
    marates.tr14(i)=tr14.Coefficients.Estimate(2)*10;
    marates.tr15(i)=tr15.Coefficients.Estimate(2)*10;


    xlim([min(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=1995)))-0.1,max(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye<=2022)))+0.1])

    legend show

    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=20; ax.FontWeight='Bold';   ax.FontName='Times New Roman'; %grid minor; grid on;
    set(gca,'fontname','Times New Roman','FontSize',20);

    ylabel('DT[cm]','FontSize',24,'FontWeight','bold');

    leg=legend(h,'location','SouthEast');
    title(leg,TGinfo.TG(i))

    xticks(1996:2:2022)
    xlim([1994 2023.5])
    pbaspect([1 .25 .25])

end

%% 

% figure 6
close all
marates2=table();
for i=1:length(tg)

    figure(i)
    % 13 month moving avarage
    T1 = length(sa_mean_m{i}.dtsa);
    sW13 = [1/24; repmat(1/12,11,1); 1/24];
    yS = conv(sa_mean_m{i}.dtsa,sW13,'same');
    yS(1:6) = yS(7); yS(T1-5:T1) = yS(T1-6);

    T2 = length(tg_mean_m{i}.dttg);
    yT = conv(tg_mean_m{i}.dttg,sW13,'same');
    yT(1:6) = yT(7); yT(T2-5:T2) = yT(T2-6);

    trtg=fitlm(decyear(tg_mean_m{i}.time),tg_mean_m{i}.dttg);
    trsa=fitlm(decyear(sa_mean_m{i}.time),sa_mean_m{i}.dtsa);


    h(1)=plot(decyear(sa_mean_m{i}.time),yS,'-r','LineWidth',1.5,'DisplayName','SA');
    hold on
    h(2)=plot(decyear(tg_mean_m{i}.time),yT,'-b','LineWidth',1.5,'DisplayName','TG');


    % h(1)=plot(decyear(sa_mean_m{i}.time),yS,'-r','LineWidth',1.5,'DisplayName',strcat('SA trend= ',num2str(trsa.Coefficients.Estimate(2)*10,2),' mm/year'));
    % hold on
    % h(2)=plot(decyear(tg_mean_m{i}.time),yT,'-b','LineWidth',1.5,'DisplayName',strcat('TG trend=',num2str(trtg.Coefficients.Estimate(2)*10,2),' mm/year'));

    % plot(decyear(tg_mean_m{i}.time),trtg.Fitted,'-','color',[0 0.4470 0.7410],'LineWidth',1)
    % plot(decyear(sa_mean_m{i}.time),trsa.Fitted,'-','color',[0.6350 0.0780 0.1840],'LineWidth',1)



    tr1=fitlm(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2000&tg_mean_m{i}.ye<2010)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=2000&tg_mean_m{i}.ye<2010));
    tr3=fitlm(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2010&tg_mean_m{i}.ye<2020)),tg_mean_m{i}.dttg(tg_mean_m{i}.ye>=2010&tg_mean_m{i}.ye<2020));

    tr11=fitlm(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=2000&sa_mean_m{i}.ye<2010)),sa_mean_m{i}.dtsa(sa_mean_m{i}.ye>=2000&sa_mean_m{i}.ye<2010));
    tr13=fitlm(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=2010&sa_mean_m{i}.ye<2020)),sa_mean_m{i}.dtsa(sa_mean_m{i}.ye>=2010&sa_mean_m{i}.ye<2020));


    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2000&tg_mean_m{i}.ye<2010)),tr1.Fitted,'-','color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
    plot(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=2010&tg_mean_m{i}.ye<2020)),tr3.Fitted,'-','color',[0.4940 0.1840 0.5560],'LineWidth',1.5)


    plot(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=2000&sa_mean_m{i}.ye<2010)),tr11.Fitted,'--','color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
    plot(decyear(sa_mean_m{i}.time(sa_mean_m{i}.ye>=2010&sa_mean_m{i}.ye<2020)),tr13.Fitted,'--','color',[0.4940 0.1840 0.5560],'LineWidth',1.5)

    % TG
    marates2.tr1(i)=tr1.Coefficients.Estimate(2)*10;
    marates2.tr3(i)=tr3.Coefficients.Estimate(2)*10;
    % SA
    marates2.tr11(i)=tr11.Coefficients.Estimate(2)*10;
    marates2.tr13(i)=tr13.Coefficients.Estimate(2)*10;

    xlim([min(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye>=1995)))-0.1,max(decyear(tg_mean_m{i}.time(tg_mean_m{i}.ye<=2022)))+0.1])

    legend show

    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=20; ax.FontWeight='Bold';   ax.FontName='Times New Roman'; %grid minor; grid on;
    set(gca,'fontname','Times New Roman','FontSize',20);

    ylabel('DT[cm]','FontSize',24,'FontWeight','bold');

    leg=legend(h,'location','SouthEast');
    title(leg,TGinfo.TG(i))

    xticks(1996:2:2022)
    xlim([1994 2023.5])
    pbaspect([1 .25 .25])

end

%% decadal TG trend line plot
% figure x3
close all
limit=5;


for i=1:length(tg)
    
    figure(i)
    
    mintime=min(sa{i}.t);
    maxtime=max(sa{i}.t);
    
    ti1=decyear(sa{i}.t(sa{i}.dist>limit&sa{i}.year_sa<2010&sa{i}.year_sa>=2000));
    ti2=decyear(sa{i}.t(sa{i}.dist>limit&sa{i}.year_sa<2020&sa{i}.year_sa>=2010));
    ti4=decyear(sa{i}.t(sa{i}.dist>limit));
    
    dtsa1=sa{i}.dt2(sa{i}.dist>limit&sa{i}.year_sa<2010&sa{i}.year_sa>=2000);
    dtsa2=sa{i}.dt2(sa{i}.dist>limit&sa{i}.year_sa<2020&sa{i}.year_sa>=2010);
    dtsa4=sa{i}.dt2(sa{i}.dist>limit);
    
    %TG
    dttg1=tg{i}.dt3(tg{i}.y<2010&tg{i}.y>=2000);
    dttg2=tg{i}.dt3(tg{i}.y<2020&tg{i}.y>=2010);
    dttg4=tg{i}.dt3(tg{i}.Time<=maxtime&tg{i}.Time>=mintime);
    
    tgti1=decyear(tg{i}.Time(tg{i}.y<2010&tg{i}.y>=2000));
    tgti2=decyear(tg{i}.Time(tg{i}.y<2020&tg{i}.y>=2010));
    tgti4=decyear(tg{i}.Time(tg{i}.Time<=maxtime&tg{i}.Time>=mintime));
    
    %trend
    tr_tg1=fitlm(tgti1,dttg1);
    tr_tg2=fitlm(tgti2,dttg2);
    tr_tg4=fitlm(tgti4,dttg4);
    
    tr_sa1=fitlm(ti1,dtsa1);
    tr_sa2=fitlm(ti2,dtsa2);
    tr_sa4=fitlm(ti4,dtsa4);
    
    
    h(1)=plot(tgti4,dttg4,'.','color',[.7 .7 .7],'DisplayName',strcat('ASL_T_G, trend_{1995-2022}: ',num2str(tr_tg4.Coefficients.Estimate(2)*10,2),' mm/year'));
    hold on
    h(2)=plot(tgti1,tr_tg1.Fitted,'--b','LineWidth',1.5,'DisplayName',strcat('TG trend_{2000-2010}: ',num2str(tr_tg1.Coefficients.Estimate(2)*10,2),' mm/year'));
    h(3)=plot(tgti2,tr_tg2.Fitted,'-b','LineWidth',1.5,'DisplayName',strcat('TG trend_{2010-2020}: ',num2str(tr_tg2.Coefficients.Estimate(2)*10,2),' mm/year'));
   
    h(4)=plot(tgti4,tr_tg4.Fitted,'-r','LineWidth',1.5,'DisplayName',strcat('TG trend_{1995-2022}: ',num2str(tr_tg4.Coefficients.Estimate(2)*10,2),' mm/year'));

    
    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=20; ax.FontWeight='Bold';   ax.FontName='Times New Roman'; %grid minor; grid on;
    set(gca,'fontname','Times New Roman','FontSize',20);
    
    ylabel('DT[cm]','FontSize',24,'FontWeight','bold');
    
    % leg=legend(h(1:3),'location','NorthWest');
     leg=legend(h(1:4),'location','NorthWest');

    title(leg,TGinfo.TG(i))
    
    xticks(1996:2:2022)
    xlim([1994 2023.5])
    pbaspect([1 .25 .25])
    
    clearvars h leg ax tr_sa4 tr_sa2 tr_sa1 tr_tg4 tr_tg2 tr_tg1 tgti4 tgti2 tgti1 dtsa1 dtsa2 dtsa4 ti1 ti2 ti4 maxtime mintime
end

%% moving avarge

ft = fittype('a1+a2*x+a*cos(pi*x)+b*sin(pi*x)+c*cos(2*pi*x)+d*sin(2*pi*x)');
limit=5;
for i=1:4


    mintime=min(sa{i}.t);
    maxtime=max(sa{i}.t);
    sa{i}.mo=month(sa{i}.t);

    %SA
    x1=decyear(sa{i}.t(sa{i}.dist>limit));
    y1=sa{i}.dt(sa{i}.dist>limit)*100;

    %TG
    x2=decyear(tg{i}.Time(tg{i}.Time<=maxtime&tg{i}.Time>=mintime));
    y2=tg{i}.dt3(tg{i}.Time<=maxtime&tg{i}.Time>=mintime);

    fitted1=fit(x1,y1,ft);
    % fitted2=fit(x2,y2,ft);

    tr1=fitlm(x1,y1);
    tr2=fitlm(x2,y2);

    mo_sa=movmean(sa_mean_m{i}.dtsa*100,12,'SamplePoints',sa_mean_m{i}.tm);
    mo_tg=movmean(tg_mean_m{i}.dttg,12,'SamplePoints',tg_mean_m{i}.tm);

    % sm_sa=smoothdata(sa_mean_m{i}.dtsa*100,'movmean',120,'SamplePoints',sa_mean_m{i}.tm);


    [lt_sa,st_sa,r_sa] = trenddecomp(sa_mean_m{i}.dtsa*100);
    lt_tg=trenddecomp(tg_mean_m{i}.dttg);


    ac_sa = tsaccel(lt_sa);
    ac_tg = tsaccel(lt_tg);

    acsa=fitlm(x1,tr1.Fitted);
    actg=fitlm(x2,tr2.Fitted);
    plot(x1,acsa.Fitted)
    hold on
    plot(x2,actg.Fitted)

    plot(decyear(sa_mean_m{i}.time),ac_sa)
    hold on
    plot(decyear(tg_mean_m{i}.time),ac_tg)

    actrsa=fitlm(decyear(sa_mean_m{i}.time),ac_sa);
    actrtg=fitlm(decyear(tg_mean_m{i}.time),ac_tg);

    plot(decyear(sa_mean_m{i}.time),actrsa.Fitted)
    hold on
    plot(decyear(tg_mean_m{i}.time),actrtg.Fitted)




    sst1=fitted1.a1+fitted1.a2*x1+fitted1.a*cos(pi*x1)+fitted1.b*sin(pi*x1)+fitted1.c*cos(2*pi*x1)+fitted1.d*sin(2*pi*x1);
    % sst2=fitted2.a1+fitted2.a2*x2+fitted2.a*cos(pi*x2)+fitted2.b*sin(pi*x2)+fitted2.c*cos(2*pi*x2)+fitted2.d*sin(2*pi*x2);


    % 13 month moving avarage
    T1 = length(sa_mean_m{i}.dtsa*100);
    sW13 = [1/24; repmat(1/12,11,1); 1/24];
    yS = conv(sa_mean_m{i}.dtsa*100,sW13,'same');
    yS(1:6) = yS(7); yS(T1-5:T1) = yS(T1-6);

    % % find seasonality
    % Fs = 12;
    % T = 1/Fs;
    % t = (0:T1-1)*T;

    figure(i)
    %SA

    plot(decyear(sa_mean_m{i}.time),mo_sa,'-r','LineWidth',1.5)
    hold on
    h(1)=plot(x1,tr1.Fitted,'-','color',[0.6350 0.0780 0.1840],'DisplayName',strcat('SA_{Linear Trend}: ',num2str(tr1.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',3);
    plot(decyear(sa_mean_m{i}.time),lt_sa,'-','color',[.5 .5 .5],'LineWidth',3,'DisplayName','long trend')
    plot(decyear(sa_mean_m{i}.time),yS,'-','color',[0.4660 0.6740 0.1880],'LineWidth',1.5,'DisplayName','13-month Moving Average');
    % plot(decyear(sa{i}.t(sa{i}.dist>limit)),sst1,'-','color',[0.4940 0.1840 0.5560],'LineWidth',1.5,'DisplayName','Seasonality');
    % plot(decyear(sa_mean_m{i}.time),st_sa,'-k','LineWidth',1.5,'DisplayName','Seasonality');


    %TG
    plot(decyear(tg_mean_m{i}.time),mo_tg,'-','color',[0 0.4470 0.7410],'LineWidth',1.5)
    h(2)=plot(x2,tr2.Fitted,'-b','DisplayName',strcat('TG_{Linear Trend}: ',num2str(tr2.Coefficients.Estimate(2)*10,2),' mm/year'),'LineWidth',3);
    plot(decyear(tg_mean_m{i}.time),lt_tg,'--','color',[.5 .5 .5],'LineWidth',3)


    ylabel('DT[cm]','FontSize',24,'FontWeight','bold');

    % plot(x1,y1,'.','color',[.7 .7 .7])
    % plot(x1,sst1,'-r','linewidth',1.5)
    % plot(x2,sst2,'-b','linewidth',1.5)
    leg=legend(h,'location','SouthEast');
    title(leg,TGinfo.TG(i))

    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; ax.FontName='Times New Roman';
    xticks(1996:2:2022)
    xlim([1994 2023.5])
    pbaspect([1 .25 .25])

    clearvars ax leg h sst1 lt_sa mo_tg mo_sa   tr1 tr2 fitted1 x1 x2 y1 y2 mintime maxtime
end
%% EMD decompostion
for k=1:13
decom{k}=table();

[m,~]=find(isnan(mean_all{k}.sa));
mean_all{k}(m,:)=[];
[m,~]=find(isnan(mean_all{k}.tg));
mean_all{k}(m,:)=[];

[imf1,residual1,~] = emd(mean_all{k}.sa,'MaxNumIMF',5);
[imf2,residual2,~] = emd(mean_all{k}.tg,'MaxNumIMF',5);

residual1=residual1-(mean(residual1)-mean(residual2));


decom{k}.imfsa1=imf1(:,1);
decom{k}.imfsa2=imf1(:,2);
decom{k}.imfsa3=imf1(:,3);
decom{k}.imfsa4=imf1(:,4);
decom{k}.imfsa5=imf1(:,5);
decom{k}.ressa=residual1;

decom{k}.imftg1=imf2(:,1);
decom{k}.imftg2=imf2(:,2);
decom{k}.imftg3=imf2(:,3);
decom{k}.imftg4=imf2(:,4);
decom{k}.imftg5=imf2(:,5);
decom{k}.restg=residual2;

decom{k}.time=mean_all{k}.time;

clearvars residual1 residual2 imf1 imf2 m

end


%% plot EMD decomposition
for i=1:6
subplot(6,1,i)
if i~=6
plot(mean_all{k}.time,imf2(:,i),'LineWidth',1.5)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;
grid minor
ylabel(strcat('IMF',num2str(i)))
xticks(datetime('01-Jan-1995'):calmonths(36):datetime('01-Jan-2023'))
xticklabels([])

else
plot(mean_all{k}.time,residual2,'LineWidth',1.5)
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on; 
xticks(datetime('01-Jan-1995'):calmonths(36):datetime('01-Jan-2023'))
ylabel('Residual')
grid minor
end

end



%% EMD plot
% figure 7

close all
for k=1:13
figure(k)


dttg=tg{k}.dt3(tg{k}.date<decyear(max(decom{k}.time)));
tgti=decyear(tg{i}.date(tg{k}.date<decyear(max(decom{k}.time))));
tr_tg=fitlm(tgti,dttg);




r1=corrcoef(decom{k}.restg,decom{k}.ressa);
r=r1(1,2)*100;

plot(decyear(decom{k}.time),decom{k}.ressa,'-r','LineWidth',1.5,'DisplayName',strcat('SA . ','RMSE=',num2str(rms(decom{k}.ressa-decom{k}.restg),3),'cm'))
hold on
plot(decyear(decom{k}.time),decom{k}.restg,'-b','LineWidth',1.5,'DisplayName',strcat('TG . ','R=',num2str(r,2)))
plot(tgti,tr_tg.Fitted,'--k','LineWidth',1.5,'DisplayName','Linear Trend')



leg = legend('show','location','southeast');
title(leg,TGinfo.TG(k))

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; ax.FontName='Times New Roman'; grid minor; grid on
xticks(1996:2:2022)
xlim([1994 2023.5])
pbaspect([1 .25 .25])

end
%% statistics of EMD
% table 4

for k=1:13

 
r1=corrcoef(decom{k}.restg,decom{k}.ressa);
r=r1(1,2)*100;
recomresult(k,1)=rms(decom{k}.ressa-decom{k}.restg);
recomresult(k,2)=r;


    
    
 tr1=fitlm(decyear(decom{k}.time),decom{k}.ressa);
 recomresult(k,3)=tr1.Coefficients.Estimate(2)*10;
 
     
    
 tr2=fitlm(decyear(decom{k}.time),decom{k}.restg);
 recomresult(k,4)=tr2.Coefficients.Estimate(2)*10;

end

%% trend from the coast 
limit=5;

% The 100(1 – α)% confidence intervals for regression coefficient
alpha=0.05;
k=1;

for i=1:length(tg)
    rescoast{i}=table();

    while limit<55
    %SA
    ti=decyear(sa{i}.t(sa{i}.dist>5&sa{i}.dist<=limit+5));
    dtsa=sa{i}.dt2(sa{i}.dist>5&sa{i}.dist<=limit+5);
    tr_sa=fitlm(ti,dtsa);
    rescoast{i}.dtsa(k,1)=tr_sa.Coefficients.Estimate(2)*10;
    %  coefficient standard errors
    co_sa=coefCI(tr_sa,alpha);
    rescoast{i}.ersa(k,1)=(co_sa(2,2)-co_sa(2,1))/2;
    rescoast{i}.dis(k,1)=limit;
    rescoast{i}.siz(k,1)=size(dtsa,1);
    limit=limit+5;
    k=k+1;
    clearvars co_sa ti dtsa tr_sa
    end
    k=1; limit=5;
end

clear k


% 
% rescoast{1, 1}.dtsa(9) = nan;
% rescoast{1, 1}.dtsa(10) = nan;
% rescoast{1, 11}.dtsa(10) = nan;
% rescoast{1, 12}.dtsa(9) = nan;
% rescoast{1, 12}.dtsa(10) = nan;
% rescoast{1, 9}.dtsa(10) = nan;
% rescoast{1, 9}.ersa(9) = 0.078011547516259100;
% rescoast{1, 10}.dtsa(10) = 7.063740733607847;
% rescoast{1, 11}.ersa(9) = 0.0321945645786299;
% rescoast{1, 12}.ersa(8) = 0.0186691520830795;
% rescoast{1, 1}.ersa(8) = 0.0178853498454715;
% rescoast{1, 6}.ersa(10) = 0.0288054680804690;
% rescoast{1, 6}.ersa(9) = 0.0186047798418315;
% rescoast{1, 8}.dtsa(8) = -4.780453869958220;
% rescoast{1, 8}.dtsa(10) = 8.127959683111893;
% rescoast{1, 8}.ersa(8) = 0.1149200520154485;
% rescoast{1, 8}.ersa(9) = 0.3614457440822240;
% rescoast{1, 9}.dtsa(9) = nan;
% rescoast{1, 9}.ersa(9) = nan;
% rescoast{1, 8}.ersa(10) = 0.0964888340389195;

%% plot trend from the coast
% Figure 11
close all


for k=1:13
    figure(k)
    
    
    
    yyaxis right
    bar(rescoast{k}.dis,rescoast{k}.siz/1000,0.2,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',0.1)
    %  ylim([0 2000000])
    %  yticks(1:200000:600000)
    %  plot(rescoast{k}.dis,rescoast{k}.siz,'s-k')
    ylabel('No. data')
    ax=gca;
    ax.YAxis(2).Color = [0.4940 0.1840 0.5560];
    
    ytickformat('%g T')
    
    hold on
    
    yyaxis left
    plot(rescoast{k}.dis,rescoast{k}.dtsa,'--r')
    hold on
    h(1)=plot(rescoast{k}.dis,rescoast{k}.dtsa,'or','MarkerFaceColor','r','MarkerSize',10,'DisplayName','SA');
    
    for i=1:10
        
        plot(rescoast{k}.dis(i),[rescoast{k}.dtsa(i)-3*(rescoast{k}.ersa(i)*10);rescoast{k}.dtsa(i)+3*((rescoast{k}.ersa(i)*10))],'_-r','LineWidth',3)
        
        plot([rescoast{k}.dis(i);rescoast{k}.dis(i)],[rescoast{k}.dtsa(i)-3*(rescoast{k}.ersa(i)*10);rescoast{k}.dtsa(i)+3*(rescoast{k}.ersa(i)*10)],'-r','LineWidth',3)
    end
    
    plot(0,[res.dttg4(k)-3*(res.ertg4(k)*10);res.dttg4(k)+3*(res.ertg4(k)*10)],'_-b','LineWidth',3)
    plot([0;0],[res.dttg4(k)-3*(res.ertg4(k)*10);res.dttg4(k)+3*(res.ertg4(k)*10)],'-b','LineWidth',3)
    h(2)=plot(0,res.dttg4(k),'sb','MarkerFaceColor','b','MarkerSize',10,'DisplayName','TG');
    
    ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; grid on
    set(gca,'fontname','Times New Roman');
    leg=legend(h,'location','SouthEast','NumColumns',2);
    title(leg,TGinfo.TG(k))
    
    xlim([-5 55])
    xticks(0:5:50)
    % ylim([-1 6])
    pbaspect([1 .35 .35])
    
    ylabel('Sea Level Trend [mm/year]','FontSize',24,'FontWeight','bold');
    xlabel('Distance to coast [km]','FontSize',24,'FontWeight','bold');
    ax.YAxis(1).Color = 'k';
    
 
 end

clearvars k ax leg i h 


%% 
for i=1:10
    for k=1:13
delta(k,1)=rescoast{k}.dtsa(i)-res.dttg4(k);
rmsecoast(i,1)=rms(delta);
clear delta
    end
end

%%

% Specify the parameters of a signal with a sampling frequency of 12 Hz (monthly data) and a signal duration of observations.

% Compute the Fourier transform of the signal.
a=tg_mean{1, 1}.tg_detrend;
a=a-mean(a);

Y = fft(a.*hann(numel(a)));

% Compute the single-sided amplitude spectrum of the signal.

f = Fs*(0:(L-1)/2)/L;
P2 = abs(Y/L);
P1 = P2(1:(L+1)/2);
P1(2:end) = 2*P1(2:end);

% plot the single-sided spectrum.
subplot (2,3,4:6)
xline(1,'--k','LineWidth',1.5)
xline(2,'--k','LineWidth',1.5)
hold on
plot(f,P1','color',[0.4940 0.1840 0.5560],'LineWidth',2)

[~,m]=max(P1);
plot(f(m),P1(m),'v','MarkerFaceColor','red','MarkerEdgeColor','k','MarkerSize',10)
[~,n]=find(f==2);
plot(f(n),P1(n),'v','MarkerFaceColor','red','MarkerEdgeColor','k','MarkerSize',10)




%monthly mean boxchart
boxchart(tg_mean{1, 1}.mo,tg_mean{1, 1}.tg_detrend,'MarkerStyle','none','WhiskerLineColor','none')
hold on
meanWeight = groupsummary(tg_mean{1, 1}.tg_detrend,tg_mean{1, 1}.mo,'mean');
h(1)=plot(meanWeight,'-r','LineWidth',2,'DisplayName','Monthly Mean DT_{detrend}');
% boxplot(tg_mean{1, 1}.tg_detrend,tg_mean{1, 1}.mo)
xticks(1:1:12)
xlim([0 13])
ylim([-50 100]) % check

xlabel("Month number")
legend(h(1))
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=20; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',20);
box on


title("Single-Sided Spectrum of DT_T_G")
xlabel("Period [year]")
ylabel("Power")
% xticks(0.25:.25:4)
xticks([.25 .5 1 2 2.5 4])
xticklabels(1./[.25 .5 1 2 2.5 4])

ax=gca; ax.GridAlpha = 0.3; ax.FontSize=20; ax.FontWeight='Bold'; grid on;  ax.FontName='Times New Roman'; %grid minor;
set(gca,'fontname','Times New Roman','FontSize',20);
box on
% boxplot(tg_mean{1, 1}.tg_detrend,tg_mean{1, 1}.mo)

clearvars a meanWeight h




[p1,s1,mu] = polyfit(ti,dtsa,1);
[y_fit1,f1]= polyval(p1,ti,s1);

% 95% prediction interval y±2Δ
plot(ti,y_fit1+2*f1,'--','color',[0.4940 0.1840 0.5560])
hold on
plot(ti,y_fit1-2*f1,'--','color',[0.4940 0.1840 0.5560])
plot(ti,y_fit1,'-r')
% predict
mdl = fitlm(ti,dtsa);
xci = linspace(1995-10,2022+10)';
[ypred,yci] = predict(mdl,xci,'simul',true);

plot(mdl)
hold on
plot(xci,yci,'g-')
