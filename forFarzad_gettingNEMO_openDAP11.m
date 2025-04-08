clc
% clear

USERNAME = "abarzandeh"; % Insert your username
PASSWORD = "213aAa213"; % Insert your password
server = "@my.cmems-du.eu"; % Copernicus Marine server
datasetID = "cmems_mod_bal_phy_my_P1D-m";
file_nc = strcat('https://',USERNAME,':', PASSWORD,'@my.cmems-du.eu/thredds/dodsC/cmems_mod_bal_phy_my_P1D-m');
ncdisp(file_nc)
%;return
    time=ncread(file_nc,'time');
    lon=ncread(file_nc,'lon');
    lat=ncread(file_nc,'lat');
    time_num=datenum(1900,1,1+double(time));
    time_date=datetime(1900,1,1,0,0,0)+days(time);
    nx=length(lon);
    ny=length(lat);
    x1=find(lon >14, 1, 'first');
    y1=find(lat  >54, 1, 'first');
    x2= find(lon >30.2, 1, 'first')-x1 +1;
    y2= find(lat  >65.8, 1, 'first')-y1 +1;

    mkdir FarzadNEMO2
    cd FarzadNEMO2
    
%%
for yr=1995:2022
        ti=find(time_num>datenum(yr,01,01,11,0,0),1,'first');
        tf=find(time_num>datenum(yr,12,31,11,0,0),1,'first');        
        st_3=[ x1 y1 ti];
        ed_3=[ x2 y2 tf-ti+1];
        st_4=[ x1 y1 1 ti];
        ed_4=[ x2 y2 1 tf-ti+1];

        eval(['sla',num2str(yr),'=squeeze(ncread(file_nc,',char(39),'sla',char(39),',st_3,ed_3));'])
        eval(['save sla',num2str(yr),'.mat sla',num2str(yr)])
        eval(['sla',num2str(yr),'=[];'])
        disp(['Data extraction precess for year = ',num2str(yr),' has been done!'])
end
   lon=lon(x1:x1+x2-1);
   lat=lat(y1:y1+y2-1);
   
    save lon.mat lon
    save lat.mat lat
  