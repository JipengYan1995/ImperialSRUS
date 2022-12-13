
%% Super-Localize ALL EVENTS
tic
SL_events =cell(no_frames,1);
parfor frame_k = 1:no_frames
    if FrameIdx(frame_k)==1
        disp(['Processing Frame#' num2str(frame_k)]);
        SL_events{frame_k}=LocalisationInEachFrame(DATA_MBonly(:,:,frame_k),X_interp_map,Z_interp_map,PSF_cor,Xcoef_threshold);
    end
end
elapsedTime = toc;
disp(['Super-localisation takes ' num2str(elapsedTime) 's']);


