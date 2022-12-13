

%%
black_VelCmap = zeros(256,3);
black_VelCmap(1:64,2) = linspace(1,0,64); %green
black_VelCmap(1:64,3) = ones(64,1); %blue
black_VelCmap(65:128,3) = linspace(1,0,64); %blue
black_VelCmap(129:192,1) = linspace(0,1,64); %red
black_VelCmap(193:256,1) = ones(1,64);
black_VelCmap(193:256,2) = linspace(0,1,64); %red
%%
MB_events_vx_cell = zeros(length(z_axis_super), length(x_axis_super));
MB_events_vy_cell=zeros(length(z_axis_super), length(x_axis_super));
count_matrix=zeros(length(z_axis_super), length(x_axis_super));
count_plus=zeros(length(z_axis_super), length(x_axis_super));
count_minus=zeros(length(z_axis_super), length(x_axis_super));
for  n=1:length(SL_events_linked_filtered)
    delta_z = SL_events_linked_filtered{n}(1).centroid_z - SL_events_linked_filtered{n}(end).centroid_z;
        delta_x = SL_events_linked_filtered{n}(1).centroid_x - SL_events_linked_filtered{n}(end).centroid_x;
     
    if (delta_z*z_super_res)^2+(delta_x*x_super_res)^2>=(1e-6*track_length_dis^2)
  
    for m=1:length(SL_events_linked_filtered{n})-1
        delta_z = SL_events_linked_filtered{n}(m).centroid_z - SL_events_linked_filtered{n}(m+1).centroid_z;
        delta_x = SL_events_linked_filtered{n}(m).centroid_x - SL_events_linked_filtered{n}(m+1).centroid_x;
        MB_travel_dist = sqrt((delta_z*z_super_res)^2 + (delta_x*x_super_res)^2);
        if MB_travel_dist>=vel_lim_min/frame_rate
            MB_travel_pixels = sqrt((delta_z)^2 + (delta_x)^2);
            current_centroid_z=SL_events_linked_filtered{n}(m).centroid_z;
            current_centroid_x=SL_events_linked_filtered{n}(m).centroid_x;
            if line_interp==1
                %         for d = 0:1/round(MB_travel_pixels*SL_events_linked_filtered(n).similarity):1   % plotting loop
                for d = 0:1/round(MB_travel_pixels):1   % plotting loop
                    %             MB_events_velocity_cell(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x)) = MB_events_velocity_cell(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))+sign(delta_z) *MB_travel_dist * frame_rate;%*SL_events_linked_filtered(n).similarity;%sign(delta_z) *
                    MB_events_vx_cell(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x)) =...
                        MB_events_vx_cell(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))-delta_x*x_super_res*frame_rate;%sign(delta_z) * *SL_events_linked_filtered(n).similarity
                    MB_events_vy_cell(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))= ...
                        MB_events_vy_cell(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))+delta_z*z_super_res*frame_rate;
                    count_matrix(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))=count_matrix(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))+1;%SL_events_linked_filtered(n).similarity;
                    if sign(delta_z)>0
                        count_plus(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))=count_plus(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))+1;%
                    elseif sign(delta_z)<0
                        count_minus(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))=count_minus(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))+1;%
                    else
                        count_minus(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))=count_minus(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))+1;%
                        count_plus(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))=count_plus(round(current_centroid_z - d*delta_z), round(current_centroid_x - d*delta_x))+1;%
                    end
                end
            else
                MB_events_vx_cell((current_centroid_z), round(current_centroid_x)) = MB_events_vx_cell((current_centroid_z), round(current_centroid_x))-delta_x*x_super_res*frame_rate;
                count_matrix((current_centroid_z), round(current_centroid_x))=count_matrix((current_centroid_z), (current_centroid_x))+1;
                MB_events_vy_cell(round(current_centroid_z), round(current_centroid_x))= ...
                    MB_events_vy_cell(round(current_centroid_z), round(current_centroid_x ))+delta_z*z_super_res*frame_rate;
                if sign(delta_z)>0
                    count_plus((current_centroid_z), round(current_centroid_x))=count_plus((current_centroid_z), round(current_centroid_x))+1;%
                elseif sign(delta_z)<0
                    count_minus((current_centroid_z), round(current_centroid_x))=count_minus((current_centroid_z), round(current_centroid_x))+1;%
                else
                    count_minus((current_centroid_z), round(current_centroid_x))=count_minus((current_centroid_z), round(current_centroid_x))+1;%
                    count_plus((current_centroid_z), round(current_centroid_x))=count_plus((current_centroid_z), round(current_centroid_x))+1;%
                end
            end
        end
    end
    end
end

MB_events_vx=MB_events_vx_cell;
MB_events_vy=MB_events_vy_cell;
clear MB_events_vx_cell MB_events_vy_cell
%% location map
sigma=mean([sigma1,sigma2])/sigma_factor;
SL_MB_events_smooth = imgaussfilt(count_matrix, sigma);
aa=zeros(101,101);
aa(51,:)=1;
aa=imgaussfilt(aa, sigma);
SL_MB_events_smooth=SL_MB_events_smooth./max(aa(:));

% Define a new colormap
temp = jet(256);
jet_black_bg = [zeros(4,3); temp];
if ~exist('colorbar_lim')
    colorbar_lim = round(max(SL_MB_events_smooth(:)))/2.5;
end
figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
imagesc(x_axis_super*1e3, z_axis_super*1e3, SL_MB_events_smooth, [0 colorbar_lim])
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
formatOut = 'dd-mm-yyyy_HHMM';
timestamp = datestr(now,formatOut);
filname = ([ result_path 'SR_Image_' timestamp ]);
print(gcf,filname,outputFormat,'-r1200')
%% direction map
SL_MB_events_smooth_plus = imgaussfilt(count_plus, sigma);
SL_MB_events_smooth_plus = max(count_plus(:))* SL_MB_events_smooth_plus/ max(SL_MB_events_smooth_plus(:));

SL_MB_events_smooth_minus = imgaussfilt(count_minus, sigma);
SL_MB_events_smooth_minus = -max(count_minus(:))* SL_MB_events_smooth_minus/ max(SL_MB_events_smooth_minus(:));
sign_val=SL_MB_events_smooth_plus+SL_MB_events_smooth_minus;
SL_MB_events_smooth_combine=sign_val;%max(SL_MB_events_smooth_plus,(SL_MB_events_smooth_minus)).*(sign(SL_MB_events_smooth_plus+SL_MB_events_smooth_minus));


clim=max(abs(SL_MB_events_smooth_combine(:)))/3;

figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
imagesc(x_axis_super*1e3, z_axis_super*1e3,SL_MB_events_smooth_combine,[-clim clim])
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
filname = ([ result_path 'SR_Image_direction' timestamp]);
print(gcf,filname,outputFormat,'-r1200')
%% seperated direction map 

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
filname = ([ result_path 'SR_Image_direction_up' timestamp]);
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
filname = ([ result_path 'SR_Image_direction_down' timestamp]);
print(gcf,filname,outputFormat,'-r1200')
%% velocitymap
sigma_vel=sigma;
%  %% circle

f_fun=@(x) (exp(-x^2/2/sigma^2)-0.05);
radius=round(fzero(f_fun,[0 sigma*20]));

conv_kernel=strel('disk',radius);
conv_kernel=conv_kernel.Neighborhood;
MB_events_velocity=sqrt(MB_events_vy.^2+MB_events_vx.^2);
MB_events_velocity_log=count_matrix;
MB_events_velocity_log_smooth=conv2(MB_events_velocity_log,conv_kernel,'same');

MB_events_velocity_smooth_vx=conv2(MB_events_vx,conv_kernel,'same')./MB_events_velocity_log_smooth;
MB_events_velocity_smooth_vx(isnan(MB_events_velocity_smooth_vx))=0;
MB_events_velocity_smooth_vy=conv2(MB_events_vy,conv_kernel,'same')./MB_events_velocity_log_smooth;
MB_events_velocity_smooth_vy(isnan(MB_events_velocity_smooth_vy))=0;
MB_events_velocity_smooth=sqrt(MB_events_velocity_smooth_vx.^2+MB_events_velocity_smooth_vy.^2);
MB_events_velocity_draw=MB_events_velocity_smooth;%sqrt(MB_events_velocity_smooth_vx.^2+MB_events_velocity_smooth_vy.^2);

MB_events_angle= atan2(MB_events_velocity_smooth_vy,MB_events_velocity_smooth_vx);
MB_events_angle=MB_events_angle/pi*180;
MB_events_angle_smooth=MB_events_angle;

%% absolute velocity with brightness

if ~exist('colorbar_lim_vel')
    colorbar_lim_vel=max(MB_events_velocity_smooth(:))/1.5;
end

sat=linspace(1,1,256);
light=linspace(0,1,256);

H=[linspace(240,0,256)/360];

HSVmap=ones(256,256,3);
HSVmap(:,:,1)=repmat(H',1,256);
HSVmap(:,:,3)=repmat(light,256,1);
HSVmap(:,:,2)=repmat([sat]',1,256);

RGBmapABS=hsv2rgb(HSVmap);
colorbar_lim_draw=colorbar_lim/2;

% %%
MB_events_vy_dis=conv2(MB_events_vy,conv_kernel,'same');
MB_events_velocity_smooth_lim=MB_events_velocity_draw;
MB_events_velocity_smooth_lim(MB_events_velocity_smooth_lim>colorbar_lim_vel)=colorbar_lim_vel;
MB_events_velocity_smooth_lim([isnan(MB_events_velocity_smooth_lim)])=0;
AbsVelMapCData=zeros(size(MB_events_velocity_smooth_lim,1),size(MB_events_velocity_smooth_lim,2),3);
MB_events_velocity_smooth_recale=round(MB_events_velocity_smooth_lim/colorbar_lim_vel*255)+1;%round(rescale(MB_events_velocity_smooth_lim,1,256));
MB_events_localization_rescale=SL_MB_events_smooth;
MB_events_localization_rescale(MB_events_localization_rescale>colorbar_lim_draw)=colorbar_lim_draw;
MB_events_localization_rescale=round(MB_events_localization_rescale/colorbar_lim_draw*255)+1;%round(rescale(MB_events_localization_rescale,1,256));
% %%
for i=1:size(AbsVelMapCData,1)
    for j=1:size(AbsVelMapCData,2)
          AbsVelMapCData(i,j,:)=RGBmapABS(MB_events_velocity_smooth_recale(i,j),MB_events_localization_rescale(i,j),:);
    end
end
% %%
figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
subplot('Position',[0.1 0.1 0.7 0.8])
image('XData',(x_axis_super*1e3), 'YData',(z_axis_super*1e3),'CData',(AbsVelMapCData))
axis image
set(gca,'YDir','reverse')
xlabel('Lateral (mm)')
ylabel('Depth (mm)')
title(['Velocity Map '])
set(gca,'FontSize',14,'Fontname','Arial')



subplot('Position',[0.85 0.1 0.1 0.8])
imagesc('XData',(0:255)/255*colorbar_lim_draw,'YData',(0:255)/255*colorbar_lim_vel*1e3,'CData',RGBmapABS)
xlabel('MB density')
ylabel('Velocity (mm/s)')
title('Velocity colormap')
axis tight
set(gca,'FontSize',8,'Fontname','Arial')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
timestamp = datestr(now,formatOut);
filname = ([result_path 'SR_vel_absolute_brightness' timestamp ]);
print(gcf,filname,outputFormat,'-r1200')
%% Angle map with amplitude
sat=linspace(1,1,256);
light=linspace(0,1,256);
H=[linspace(90,360,270)/360 linspace(0,90,91)/360 ];
HSVmap=ones(361,256,3);
HSVmap(:,:,1)=repmat(H',1,256);
HSVmap(:,:,3)=repmat(light,361,1);
HSVmap(:,:,2)=1;
RGBmapAngle=hsv2rgb(HSVmap);
%%
Circle_colormap=ones(401,401,3);
[x_mesh,y_mesh]=meshgrid(-200:1:200,200:-1:-200);
circle_angle= atan2(y_mesh,x_mesh);
circle_angle=circle_angle/pi*180;
circle_angle=rescale(circle_angle,1,361);
for i=1:size(Circle_colormap,1)
    for j=1:size(Circle_colormap,2)
        if sqrt(x_mesh(i,j)^2+y_mesh(i,j)^2)<=200
            Circle_colormap(i,j,:)=RGBmapAngle(round(circle_angle(i,j)),round(sqrt(x_mesh(i,j)^2+y_mesh(i,j)^2)/200*255+1),:);
        end
    end
end
%%
AngleMapCData=zeros(size(MB_events_angle_smooth,1),size(MB_events_angle_smooth,2),3);
MB_events_angle_smooth_rescale=round(MB_events_angle_smooth)+181;%rescale(MB_events_angle_smooth,1,361);
for i=1:size(AngleMapCData,1)
    for j=1:size(AngleMapCData,2)
        AngleMapCData(i,j,:)=RGBmapAngle(round(MB_events_angle_smooth_rescale(i,j)),MB_events_localization_rescale(i,j),:);
    end
end
% %%
figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
subplot('Position',[0.1 0.1 0.7 0.8])
image('XData',(x_axis_super*1e3), 'YData',(z_axis_super*1e3),'CData',(AngleMapCData))
axis image
set(gca,'YDir','reverse')
xlabel('Lateral (mm)')
ylabel('Depth (mm)')
title(['Angle Map '])
set(gca,'FontSize',14,'Fontname','Arial')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
subplot('Position',[0.85 0.1 0.1 0.8])
imshow((Circle_colormap))

timestamp = datestr(now,formatOut);
filname = ([result_path 'SR_angle_brightness' timestamp ]);
print(gcf,filname,outputFormat,'-r1200')
