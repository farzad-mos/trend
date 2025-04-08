%% EMD
[imf3,residual3,~] = emd(forcast_imf{k}.tgdt,'MaxNumIMF',3);
forcast_imf{k}.imftg1=imf3(:,1);
forcast_imf{k}.imftg2=imf3(:,2);
forcast_imf{k}.imftg3=imf3(:,3);
forcast_imf{k}.restg=residual3;


for i=1:height(forcast_imf{k})
    forcast_imf{1, 1}.dgdecom(i)=imf3(i,1)+imf3(i,2)+imf3(i,3);
end


figure(2)
for i=1:size(imf,2)+1
subplot(size(imf,2)+1,1,i)

if i~=4
    plot(imf2(:,i))
    hold on
    plot(forcast_imf{1, 1}.("imf"+i),'displayname',strcat('RMSE=',num2str(rms(imf2(:,i)-forcast_imf{1, 1}.("imf"+i)))))
else
    
    plot(forcast_imf{1, 1}.tgdt)
    hold on
    plot(forcast_imf{1, 1}.sadecom,'displayname',strcat('RMSE=',num2str(rms(forcast_imf{1, 1}.dgdecom-forcast_imf{1, 1}.sadecom))))
    
end

if i~=size(imf,2)
    xticklabels([])
end

legend show
ax=gca; ax.GridAlpha = 0.3; ax.FontSize=18; ax.FontWeight='Bold';  grid on;

end
%% EWT

