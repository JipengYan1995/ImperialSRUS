MaxPSF=max(PSF_cor(:)); [r_correc,c_correc]=find(PSF_cor==MaxPSF);

FrameIntensity=squeeze(sum(sum(DATA_MBonly,2),1));
FrameIntensity=FrameIntensity.*FrameIdx;
[~,MaxInstensityFrame]=max(FrameIntensity);
fig3 = figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)-200 scr_size(4)-250]);
DATA_MBonly_frame = DATA_MBonly(:,:,MaxInstensityFrame);
hold off
imagesc(x_axis_roi*1e3, z_axis_roi*1e3,DATA_MBonly_frame);
%         colormap(sub1,gray);
colormap(gray);
colorbar;
set(gca,'FontSize',13,'Fontname','Arial')
title(['Contrast-mode Image - Frame #' num2str(MaxInstensityFrame)])
xlabel('Lateral (mm)')
ylabel('Depth (mm)')
axis image
hold on
tic
SL_events_temp=LocalisationInEachFrame(DATA_MBonly_frame,X_interp_map,Z_interp_map,PSF_cor,Xcoef_threshold);

for ii=1:length(SL_events_temp)
    plot(x_axis_super(SL_events_temp(ii).centroid_x)*1e3,z_axis_super(SL_events_temp(ii).centroid_z)*1e3,'+r')
end
elapsedTime=toc;
disp(['It takes around ' num2str(elapsedTime) 's to process one frame by one core']);