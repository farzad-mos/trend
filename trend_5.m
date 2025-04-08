load('F:\Analysis\BS\trend\data\sa\all\SA_1.mat')
load('TG.mat')


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
NAPlat=dms2degrees([52 22 53]);

for i=1:height(TGinfo)
T=eval(['TG' num2str(i)]);
[~,ia] = unique(T.date);
tg{i}=T(ia,:);

%TIDE system correction; zero to mean
tg{i}.zerotomean=repmat(0.29541*(sind(TGinfo.Lat(i)).^2-sind(NAPlat).^2)+0.00042*(sind(TGinfo.Lat(i)).^4-sind(NAPlat).^4),height(tg{i}),1).*100; %add to TG data
tg{i}.dt2=tg{i}.dt+tg{i}.zerotomean;

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
me=table();
[x,tm] = findgroups(sacut{k}.ym);

me.dtsa=splitapply(@mean, sacut{k}.sadt,x);
me.dttg=splitapply(@mean, sacut{k}.tgdt,x);
me.t=tm;
sa_mean{k}=me;
clearvars me x tm

end
% adding back the time
for k=1:height(TGinfo)
% a=(sa_mean{k}.t);
sa_mean{k}.ye=round((sa_mean{k}.t)/100,0);
sa_mean{k}.mo=round((((sa_mean{k}.t)/100)-round((sa_mean{k}.t)/100,0))*100,0);
end
for k=1:height(TGinfo)
sa_mean{k}.time=datetime(sa_mean{k}.ye,sa_mean{k}.mo,1);
end
%Perform Differences - detrend 
for k=1:height(TGinfo)
% sa_mean{k}.sadt_Differences=nan(height(sa_mean{k}),1);
% sa_mean{k}.sadt_Differences(2:end)= diff(sa_mean{k}.sadt_Linear);
% me.sadt_Linear=splitapply(@mean, sacut{k}.sadt_Linear,x);
sa_mean{k}.tg_detrend=detrend(sa_mean{k}.dttg);
sa_mean{k}.sa_detrend=detrend(sa_mean{k}.dtsa);
end


% TG data: find group by catagory
% for forcasting comparision

%sort by date & monthly category
for k=1:height(TGinfo)
tg{k} = sortrows(tg{k},'Time','ascend');
tg{k}.y = year(tg{k}.Time);
tg{k}.mo = month(tg{k}.Time);
tg{k}.da = day(tg{k}.Time);
tg{k}.ym=str2num([num2str(tg{k}.y),num2str(tg{k}.mo,'%0.2i')]);
tg{k}.ym2=str2num([num2str(tg{k}.y),num2str(tg{k}.mo,'%0.2i'),num2str(tg{k}.da,'%0.2i')]);
end

% monthly mean
for k=1:height(TGinfo)
me=table();
[x,tm] = findgroups(tg{k}.ym);
me.dttg=splitapply(@mean, tg{k}.dt,x);
me.t=tm;
tg_mean{k}=me;
clearvars me x tm
end
for k=1:height(TGinfo)
tg_mean{k}.ye=round((tg_mean{k}.t)/100,0);
tg_mean{k}.mo=round((((tg_mean{k}.t)/100)-round((tg_mean{k}.t)/100,0))*100,0);
tg_mean{k}.time=datetime(tg_mean{k}.ye,tg_mean{k}.mo,1);
tg_mean{k}.tg_detrend=detrend(tg_mean{k}.dttg);
end

% daily mean
for k=1:height(TGinfo)
me=table();
[x,tm] = findgroups(tg{k}.ym2);
me.dttg=splitapply(@mean, tg{k}.dt,x);
me.t=tm;
tg_mean_day{k}=me;
clearvars me x tm
end
for k=1:height(TGinfo)
tg_mean_day{k}.ye=round((tg_mean_day{k}.t)/10000,0);
tg_mean_day{k}.mo=round( (((tg_mean_day{k}.t)/10000)-round((tg_mean_day{k}.t)/10000,0))*100,0 );

tg_mean_day{k}.da=tg_mean_day{k}.t-str2num([num2str(tg_mean_day{k}.ye),num2str(tg_mean_day{k}.mo,'%0.2i'),repmat(num2str(00,'%0.2i'),height(tg_mean_day{k}),1)]);

tg_mean_day{k}.time=datetime(tg_mean_day{k}.ye,tg_mean_day{k}.mo,tg_mean_day{k}.da);
tg_mean_day{k}.tg_detrend=detrend(tg_mean_day{k}.dttg);
end
