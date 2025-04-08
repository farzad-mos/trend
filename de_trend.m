function [det,tre]=de_trend(x,y)
% input variables:
% x= time of y data in decyear format
% y= seasoned data ex: SSH or DT

p1 = polyfit(x,y,1);
tre = polyval(p1,x);
% tr1=fitlm(x,y);
% 
% a=tr1.Coefficients.Estimate(2); %p1(1)
 
det=y-tre;
 
end