

%% location map


figure
imagesc(x_axis_super*1e3, z_axis_super*1e3,bw_mask)
colormap(gray)
axis image
formatOut = 'dd-mm-yyyy_HHMM';
timestamp = datestr(now,formatOut);
filname = ([result_path 'Mask' timestamp ]);
% filname = (['SR_vel_magnitude_mask' timestamp ]);
print(gcf,filname,outputFormat,'-r1200')
save([result_path 'mask_for_calculation' timestamp ],'bw_mask')
count_matrix=count_matrix.*bw_mask;
sigma=mean([sigma1,sigma2])/sigma_factor;
SL_MB_events_smooth = imgaussfilt(count_matrix, sigma);

% %%
aa=zeros(101,101);
aa(51,:)=1;
aa=imgaussfilt(aa, sigma);

SL_MB_events_smooth=SL_MB_events_smooth./max(aa(:));

% Define a new colormap
temp = jet(256);
jet_black_bg = [zeros(4,3); temp];

%     %%
% colorbar_lim =  (max(SL_MB_events_smooth(:)))/2.5;
figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
imagesc(x_axis_super*1e3, z_axis_super*1e3, SL_MB_events_smooth, [0 colorbar_lim])
% colormap(jet_black_bg)   % colormap(h1, jet_black_bg)
colormap(hot)
colorbar
axis image
xlabel('Lateral (mm)')
ylabel('Depth (mm)')
title(['Super-Localized MBs'])
set(gca,'FontSize',14,'Fontname','Arial')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

filname = ([result_path 'SR_Image_mask' timestamp ]);
% filname = (['SR_vel_magnitude_mask' timestamp ]);
print(gcf,filname,outputFormat,'-r1200')

%% direction map
count_plus=count_plus.*bw_mask;
count_minus=count_minus.*bw_mask;

SL_MB_events_smooth_plus = imgaussfilt(count_plus, sigma);
SL_MB_events_smooth_plus = max(count_plus(:))* SL_MB_events_smooth_plus/ max(SL_MB_events_smooth_plus(:));

SL_MB_events_smooth_minus = imgaussfilt(count_minus, sigma);
SL_MB_events_smooth_minus = -max(count_minus(:))* SL_MB_events_smooth_minus/ max(SL_MB_events_smooth_minus(:));
sign_val=SL_MB_events_smooth_plus+SL_MB_events_smooth_minus;
% SL_MB_events_smooth_combine=max(SL_MB_events_smooth_plus,(SL_MB_events_smooth_minus)).*(sign(SL_MB_events_smooth_plus+SL_MB_events_smooth_minus));
black_VelCmap = zeros(256,3);
black_VelCmap(1:64,2) = linspace(1,0,64); %green
black_VelCmap(1:64,3) = ones(64,1); %blue
black_VelCmap(65:128,3) = linspace(1,0,64); %blue
black_VelCmap(129:192,1) = linspace(0,1,64); %red
black_VelCmap(193:256,1) = ones(1,64);
black_VelCmap(193:256,2) = linspace(0,1,64); %red
% sign_val;
% clim=max(abs(SL_MB_events_smooth_combine(:)))/10;
clim=max(abs(sign_val(:)))/5;
figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
imagesc(x_axis_super*1e3, z_axis_super*1e3,sign_val,[-clim clim])
% colormap(jet_black_bg)   % colormap(h1, jet_black_bg)
colormap(black_VelCmap)
colorbar
axis image
xlabel('Lateral (mm)')
ylabel('Depth (mm)')
title(['Super-Localized MBs with direction'])
set(gca,'FontSize',14,'Fontname','Arial')
filname = ([result_path 'SR_Image_with_direction_mask' timestamp ]);
% filname = (['SR_vel_magnitude_mask' timestamp ]);
print(gcf,filname,outputFormat,'-r1200')

%% spereted direction map with mask
clim=max(abs(SL_MB_events_smooth_combine(:)))/3;

figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
imagesc(x_axis_super*1e3, z_axis_super*1e3,SL_MB_events_smooth_plus,[-clim clim])
% colormap(jet_black_bg)   % colormap(h1, jet_black_bg)
colormap(black_VelCmap)
colorbar
axis image
xlabel('Lateral (mm)')
ylabel('Depth (mm)')
title(['Super-Localized MBs with direction'])
set(gca,'FontSize',14,'Fontname','Arial')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
filname = ([ result_path 'SR_Image_direction_up_mask' timestamp]);
print(gcf,filname,outputFormat,'-r1200')
figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
imagesc(x_axis_super*1e3, z_axis_super*1e3,SL_MB_events_smooth_minus,[-clim clim])
% colormap(jet_black_bg)   % colormap(h1, jet_black_bg)
colormap(black_VelCmap)
colorbar
axis image
xlabel('Lateral (mm)')
ylabel('Depth (mm)')
title(['Super-Localized MBs with direction'])
set(gca,'FontSize',14,'Fontname','Arial')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
filname = ([ result_path 'SR_Image_direction_down_mask' timestamp]);
print(gcf,filname,outputFormat,'-r1200')


%% velocity magnitude with mask
MB_events_velocity_smooth=MB_events_velocity_smooth.*bw_mask;
AbsVelMapCData=AbsVelMapCData.*repmat(bw_mask,1,1,3);
figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
image('XData',(x_axis_super*1e3), 'YData',(z_axis_super*1e3),'CData',AbsVelMapCData)
axis image
set(gca,'YDir','reverse')
xlabel('Lateral (mm)')
ylabel('Depth (mm)')
title(['Absolute Velocity Map '])
set(gca,'FontSize',14,'Fontname','Arial')


set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
timestamp = datestr(now,formatOut);
filname = ([ result_path 'SR_vel_magnitude_mask' timestamp ]);
print(gcf,filname,outputFormat,'-r1200')

%% angle map with mask
MB_events_angle_smooth=MB_events_angle_smooth.*bw_mask;
MB_events_angle_smooth=MB_events_angle_smooth./bw_mask;
%  load('circle_map')
figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
AngleMapCData=AngleMapCData.*repmat(bw_mask,1,1,3);
image('XData',(x_axis_super*1e3), 'YData',(z_axis_super*1e3),'CData',AngleMapCData)
axis image
set(gca,'YDir','reverse')
xlabel('Lateral (mm)')
ylabel('Depth (mm)')
title(['Angle Map '])
set(gca,'FontSize',14,'Fontname','Arial')


set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
timestamp = datestr(now,formatOut);
filname = ([ result_path 'SR_vel_angle_mask' timestamp ]);
% set(h,'alphadata',~isnan(MB_events_angle_smooth))
print(gcf,filname,outputFormat,'-r1200')
%% velocity parameter
ind_n0=find(MB_events_velocity_smooth>0);
velocity_nonzero=MB_events_velocity_smooth(ind_n0)*1e3;
mean_velocity=mean(velocity_nonzero);
std_velocity=std(velocity_nonzero);
max_velocity=max(velocity_nonzero);
figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
h1=histogram(velocity_nonzero,'Normalization','probability');
xlabel('velocity (mm/s)')
ylabel('Probability')
colorbar off
h1.NumBins = 20;
title({['Velocity distribution. ', 'Mean velocity: ',num2str(mean_velocity),'mm/s'];...
   ['Std velocity: ',num2str(std_velocity),'mm/s','; Max velocity: ', num2str(max_velocity),'mm/s']});
timestamp = datestr(now,formatOut);
filname = ([result_path 'SR_velocity_distribution' timestamp ]);
print(gcf,filname,outputFormat,'-r1200')

%% vessel distribution parameter
ind_nonzeros=find(SL_MB_events_smooth>0);
vessel_thresh=0.9;%(1/2/pi/sigma^2*10)/10;% 0.5; % vessel density diameter, the value increase and the blood vessel become thinner
SR_structure=(SL_MB_events_smooth>=vessel_thresh);

vessel_density=sum(SR_structure(:))/sum(bw_mask(:));
[n,r] = boxcount(SR_structure,'slope');
df = -diff(log(n))./diff(log(r));
plot(df,'-o')
inval = input('Steady range of DF? ','s');
inval=split(inval);
fd_number=mean(df(str2num(inval{1}):str2num(inval{2})));

figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
imagesc(x_axis_super*1e3, z_axis_super*1e3,SR_structure)
colormap(gray)
title(['vessel density=' num2str(round(1000*vessel_density)/10) '%, fractal number=' num2str(round(1000*fd_number)/1000) ' +/- ' num2str(round(1000*std(df(4:8)))/1000)]);
timestamp = datestr(now,formatOut);
filname = ([result_path 'SR_structure' timestamp ]);
print(gcf,filname,outputFormat,'-r1200')
%% diameter calculation
min_super_res=min(x_super_res,z_super_res);
SR_structure_square=imresize(SR_structure,[round(size(SR_structure,1)*z_super_res/min_super_res),round(size(SR_structure,2)*x_super_res/min_super_res)]); % square pixel.
edtImage = bwdist(~SR_structure_square);
skelimage=bwmorph(SR_structure_square,'skel',Inf);
figure
title('Skeleton of muscle')
imagesc(skelimage)
axis image
colormap(gray)
vessel_radius=edtImage(skelimage);
vessel_diameter=2*vessel_radius;
figure
h2=histogram(vessel_diameter,'Normalization','probability');
h2.NumBins=30;
char_cell=cell(length(h2.Values),1);
for i=1:length(h2.BinEdges)
    char_cell{i}=num2str(h2.BinEdges(i)*min_super_res*1e6);
end
xticklabels(char_cell)
xlabel('Diameter (um)')
%% tortuosity

multi_bubble_analysis_Kalman_v1;
figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
h1=histogram(tortuosity,'Normalization','probability');
xlabel('tortuosity')
ylabel('Probability')
title({['Tortuosity distribution. ', 'Mean tortuosity: ',num2str(mean(tortuosity))];...
   ['Std tortuosity: ',num2str(std(tortuosity)),'; Max tortuosity: ', num2str(max(tortuosity))]});
timestamp = datestr(now,formatOut);
filname = ([result_path 'tortuosity_distribution' timestamp ]);
print(gcf,filname,outputFormat,'-r1200')

%% save and prinf parameters
save([result_path timestamp 'calculated_parameters' ],'tortuosity','vessel_density','vessel_diameter','fd_number','MB_events_vy','MB_events_vx','count_matrix','bw_mask')
fprintf('Mean velocity: %f mm/s \n', mean_velocity);
fprintf('Std velocity: %f mm/s \n', std_velocity);
fprintf('Max velocity: %f mm/s \n', max_velocity);
fprintf('Final vessel density: %f \n', vessel_density(end));
fprintf('Mean tortuosity: %f  \n', mean(tortuosity));
fprintf('Std tortuosity: %f  \n', std(tortuosity));
fprintf('Max tortuosity: %f \n', max(tortuosity));
fprintf('fractal number: %f \n', (round(1000*fd_number)/1000));
% =' num2str(round(1000*fd_number)/1000)
fprintf('Mean diameter: %f um \n', mean(vessel_diameter*x_super_res*1e6));
fprintf('Std diameter: %f  um \n', std(vessel_diameter*x_super_res*1e6));
fprintf('Max diameter: %f um \n', max(vessel_diameter*x_super_res*1e6));
