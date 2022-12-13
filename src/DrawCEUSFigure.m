%% 显示最大亮度和�?�?�的图
figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
imagesc(x_axis*1e3, z_axis*1e3, max(DATA_MBonly,[],3))
colormap(gray)
title('Max instensity Image')
axis image
formatOut = 'dd-mm-yyyy_HHMM';
timestamp = datestr(now,formatOut);
filname = ([ result_path 'Max_instensity' timestamp ]);
print(gcf,filname,outputFormat,'-r1200')

figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
imagesc(x_axis*1e3, z_axis*1e3, sum(DATA_MBonly,3))
colormap(gray)
title('Stacked Image')
axis image
timestamp = datestr(now,formatOut);
filname = ([ result_path 'Stacked_image' timestamp ]);
print(gcf,filname,outputFormat,'-r1200')