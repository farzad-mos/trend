function [des,sst]=deseason(x,y,mo)
% input variables:
% x= time of y data in decyear format
% y= seasoned data ex: SSH or DT
% mo= 12 for annual seasonlized data 
% mo= 6 for annual and semiaual data
    
 if mo==6
     ft = fittype('a*cos(pi*x)+b*sin(pi*x)+c*cos(2*pi*x)+d*sin(2*pi*x)');
 elseif mo==12
     ft = fittype('a*cos(pi*x)+b*sin(pi*x)');
 else
     dsiplay('_______')
 end
 
fitted=fit(x,y,ft);
 
 if mo==6
     sst=fitted.a*cos(pi*x)+fitted.b*sin(pi*x)+fitted.c*cos(2*pi*x)+fitted.d*sin(2*pi*x);
 elseif mo==12
     sst=fitted.a*cos(pi*x)+fitted.b*sin(pi*x);
 end
 
 des=y-sst;
 
end