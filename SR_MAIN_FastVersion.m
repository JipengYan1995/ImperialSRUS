% This code provides a simplified implementation of the framework in below paper. Deconvolution in the paper is replaced by normalised cross-correlation in this code for localisation.

%  J. Yan, T. Zhang, J. Broughton-Venner, P. Huang and M. -X. Tang,
%  "Super-Resolution Ultrasound Through Sparsity-Based Deconvolution and
%  Multi-Feature Tracking," in IEEE Transactions on Medical Imaging, vol.
%  41, no. 8, pp. 1938-1947, Aug. 2022, doi: 10.1109/TMI.2022.3152396.  

% to illustrate localisation, tracking and metric calculations
% of Ultrasound Localisation Microscopy(ULM)/ Super-resolution ultrasound
% imaging. 

% Code is published under a Creative Common license for non commercial use
% (CC-BY-NC), and therefore can be used for non-commercial, personal or
% academic use as long as the aforementioned paper is correctly cited.

clear
%% Path
scr_size = get(0, 'ScreenSize');            % get screen size(for display purposes only)
addpath(genpath('D:\src\')) %Software path
outputFormat='-djpeg'; % Save format of SR images: '-dpdf' '-djpeg' '-dpng'
data_path='F:\LabData\';% Data path
file_name='registration.mat';%Data Name
result_path=[data_path 'result\'];
mkdir(result_path)
%load data
load([data_path file_name])
%%
DATA_MBonly=registered_contrastimage;% input motion-corrected contrast image to 'DATA_MBonly'
% clear registered_contrastimage
no_frames =size(DATA_MBonly,3);         % find number of frames
FrameIdx=ones(no_frames,1);
abandoned_frames=[]; %select frames to be abandoned. or abandoned_frames=[];20:30;
FrameIdx(abandoned_frames)=0;

z_res=0.000032737;  % define pixel resolution, unit mm
x_res=z_res; % if pixel is not square, resize pixel to be square.
x_axis=((1:size(DATA_MBonly,2))-1)*x_res;  % define lateral axis
z_axis=((1:size(DATA_MBonly,1))-1)*z_res;    % define depth axis

DrawCEUSFigure  % save MIP and Stacked images
%% remove background noise
noise_thresh=0.1; % Noise threshold for removing noise
background_noise=max([mean(mean(mean(abs(DATA_MBonly),3))),noise_thresh]);
DATA_MBonly=DATA_MBonly.*(DATA_MBonly>repmat(background_noise,1,1,no_frames));%;
DATA_MBonly=DATA_MBonly/max(DATA_MBonly(:));
%show image
figure('Position',scr_size)
subplot(2,2,1)
imagesc(DATA_MBonly(:,:,floor(no_frames/5)))
colormap(gray)
subplot(2,2,2)
imagesc(DATA_MBonly(:,:,floor(2*no_frames/5)))
colormap(gray)
subplot(2,2,3)
imagesc(DATA_MBonly(:,:,floor(3*no_frames/5)))
colormap(gray)
subplot(2,2,4)
imagesc(DATA_MBonly(:,:,floor(4*no_frames/5)))
colormap(gray)

%% Define generic parameters
close all
% %% Super-resolution image parameters
SR_factor=8; % Ratio orginal pixel size to SR pixel size.

x_super_res = x_res/SR_factor;%5e-6;
z_super_res = z_res/SR_factor;%5e-6;  
SR_psf_smooth = 20e-6; % localisation uncertainty of the imaging system. needed to be measured.
sigma1 = (SR_psf_smooth/z_super_res)/2.355; % For a Guassian kernel FWHM is at 2.355 sigma (Do not change sigma, please alter SR_psf_smooth if required)
sigma2 = (SR_psf_smooth/x_super_res)/2.355; % For a Guassian kernel FWHM is at 2.355 sigma (Do not change sigma, please alter SR_psf_smooth if required)

%% Get PSF by GUI.
% Define PSF estimation parameter
patch_size=61; % pixel length for a square sub-region
x_lengthpixel=25;
z_lengthpixel=19; % pixel length of PSF
noise_level_for_gaussianfit=0.1;% threshold to cut PSF.
bubble_num_PSF=10; %number of bubbles used  for PSF estimation
PSF_Estimation_GUI
% cut PSF to smaller size to improve the speed.
PSF=abs(PSF)/max(max(abs(PSF)));
boundary_cell = bwboundaries(PSF>noise_level_for_gaussianfit,8,'holes');
boundary = boundary_cell{1};
z_boundaries = min(boundary(:,1)):max(boundary(:,1));
x_boundaries = min(boundary(:,2)):max(boundary(:,2));
PSF=PSF(z_boundaries,x_boundaries); % Make patch of PSF smaller
PSF_cor=imresize(PSF,[(size(PSF,1)-1)*round(SR_factor)+1,(size(PSF,2)-1)*round(SR_factor)+1]); %interp PSF to SR map
%% Define Super-Resolution image axes
rect_roi = [1 1 size(DATA_MBonly,2)-1 size(DATA_MBonly,1)-1];
x_axis_roi = x_axis(round(rect_roi(1)+(0:rect_roi(3))));
z_axis_roi = z_axis(round(rect_roi(2)+(0:rect_roi(4))));
x_axis_super =x_axis_roi(1) + (1:(length(x_axis_roi)-1)*SR_factor+1)*x_super_res;
z_axis_super =z_axis_roi(1) + (1:(length(z_axis_roi)-1)*SR_factor+1)*z_super_res;
[X_interp_map,Z_interp_map] = meshgrid(1:1/SR_factor:length(x_axis_roi),1:1/SR_factor:length(z_axis_roi));
%% test parameter
Xcoef_threshold=0.5; % threshold for normalised cross-correlation coefficient map
TimeEstimation
%% Super-Localise ALL Events.
SL_AllEvents_CorrParallel
%% Define Tracking parameters
vel_limit=0.01; % searching window for pairing, unit m/s.

dist_limit = vel_limit / frame_rate; % normalized 
TrackPara.Vrate=1; % initial ratio for searching window
TrackPara.Variance_R=1; % localisation noise
TrackPara.Variance_Q=2.5; % model predicition noise
TrackPara.cost_limit=1; % cost limit of graph-based assignment 
TrackPara.CostType='Position'; %cost function type' Choose 'Position',
% if you cannot rightly correct the systematic shift in localisation 
SL_events_linked = link_trajectories_Multi(SL_events, dist_limit,z_super_res,x_super_res,TrackPara);

%% filter tracking result
min_Track_length=3; % discard bubbles with appeared frames less than this number
SL_events_linked_filtered=SL_events_linked;
Track_length=zeros(length(SL_events_linked),1);
for i=1:length(SL_events_linked)
    Track_length(i)=length(SL_events_linked{i});
end
figure
histogram(Track_length)
SL_events_linked_filtered=SL_events_linked_filtered(Track_length>=min_Track_length); %
%% plot SR images
line_interp=1; % if 1, link paired bubbles with line.
colorbar_lim=10;
colorbar_lim_vel=5e-3; % define colorbar
sigma_factor=1; % Images will be blurred more with smaller value.
vel_lim_min=0e-3; % minimal accepted velocity to remove static bubbles
track_length_dis=0.1; % minimal accepted distance to remove wrong tracks.

SR_Draw_v8;


%% Draw SR images with mask
close all
C_mode=max(DATA_MBonly(round(z_axis_roi/z_res)+1,round(x_axis_roi/x_res)+1,:),[],3);
Z_interp=1:1/SR_factor:size(C_mode,1);
X_interp=1:1/SR_factor:size(C_mode,2);
[X_interpmap,Z_interpmap] = meshgrid(X_interp,Z_interp);
C_interp=interp2(C_mode,X_interpmap,Z_interpmap);
C_interp=C_interp/max(C_interp(:));
dynamic_CData=repmat(C_interp,1,1,3);
[idx,idy]=find(SL_MB_events_smooth_combine~=0);
SL_MB_events_smooth_combine_overlap=SL_MB_events_smooth_combine;
SL_MB_events_smooth_combine_overlap(SL_MB_events_smooth_combine_overlap>clim)=clim;
SL_MB_events_smooth_combine_overlap(SL_MB_events_smooth_combine_overlap<-clim)=-clim;
SL_ColorIdx=round((SL_MB_events_smooth_combine_overlap+clim)/2/clim*255)+1;
for subi=1:length(idx)
    dynamic_CData(idx(subi),idy(subi),:)=black_VelCmap(SL_ColorIdx(idx(subi),idy(subi)),:);
end
figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)/2 scr_size(4)-250])
imagesc('XData',x_axis_super*1e3,'YData', z_axis_super*1e3,'CData', dynamic_CData)
set(gca,'YDir','reverse')
axis image
xlabel('Lateral (mm)')
ylabel('Depth (mm)')
title(['Super-Localized MBs'])
set(gca,'FontSize',14,'Fontname','Arial')
bw_mask=zeros(size(SL_MB_events_smooth));
%% run this cell when you want to mask some area in
tempfig=figure('Position',scr_size);
imagesc(dynamic_CData)
axis image
title('draw mask')
roi = drawpolygon;
bw_mask_current = createMask(roi);

bw_mask=bw_mask | bw_mask_current;
close(tempfig)
%% run this cell when you want to mask some area out.
dynamic_CData=repmat(C_interp,1,1,3);
[idx,idy]=find((SL_MB_events_smooth_combine.*bw_mask)~=0);
for subi=1:length(idx)
    dynamic_CData(idx(subi),idy(subi),:)=black_VelCmap(SL_ColorIdx(idx(subi),idy(subi)),:);
end

figure('Position',scr_size)
imagesc(dynamic_CData)
title('Remove an area out')
axis image
roi = drawpolygon;
bw_mask_current = createMask(roi);
bw_mask=bw_mask.*(1-bw_mask_current);

%% calcuate various parameters
parameter_calculation_Kalman_v_cut



