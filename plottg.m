% load('TGdata.mat')
load('Z:\data\raw_data\geoid_data\NKG2015_New\NKG2015zt.mat')
% load('tglist.mat')

% a=TGinfo(:,1:4);


latlim = [53 67];
lonlim = [10 31];

ax = usamap(latlim, lonlim);
hold on
geoshow('landareas.shp','FaceColor',[0.7 0.7 0.7])
setm(ax, 'FFaceColor', [1 1 1])
setm(gca,'FLineWidth',5,'Grid','on','FontSize',18,'fontweight','bold')
contourfm(imresize(nkglat,0.6),imresize(nkglon,0.6),imresize((nkg2015),0.6),20,'edgecolor','none')

geoshow('landareas.shp','FaceColor',[0.7 0.7 0.7])

c=colorbar; 
c.Label.String = 'NKG2015_{zt} geoid [m]';
caxis([15 50])
colormap((turbo))
	
% oceanColor = [.5 .7 .9];
% setm(ax, 'FFaceColor', oceanColor)

for i=1:height(a)
    geoshow(a.Lat(i),a.Lon(i),'Marker', '^','Color','k','MarkerFaceColor','k','MarkerSize',10)
    hold on
    textm(a.Lat(i),a.Lon(i)+0.3, num2str(a.TGID(i)),'FontSize',17,'Color','k','FontWeight','Bold');
    %                                      geobasemap landcover
end

setm(gca,'FLineWidth',5,'Grid','on','FontSize',24,'fontweight','bold')


ax=gca; ax.GridAlpha = 0.3; ax.FontSize=24; ax.FontWeight='Bold'; ax.FontName='Times New Roman';
c.FontSize=24;

xLoc =1.3001e+05;
yLoc =7.0880e+06;
scaleruler('Units', 'km', 'RulerStyle', 'patches', ...
'XLoc', xLoc, 'YLoc', yLoc,"FontSize",13);