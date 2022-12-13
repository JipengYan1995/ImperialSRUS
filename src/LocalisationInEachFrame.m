function SL_events_temp=LocalisationInEachFrame(DATA_MBonly_frame,X_interp_map,Z_interp_map,PSF_cor,Gaussian_thresh)
DATA_MBonly_frame_high=interp2(DATA_MBonly_frame,X_interp_map,Z_interp_map);
MaxPSF=max(PSF_cor(:)); 
[r_correc,c_correc]=find(PSF_cor==MaxPSF);
Correlated_Map = normxcorr2(PSF_cor,DATA_MBonly_frame_high);
% MAIN LOOP - Detect MB events and discard the false events
Correlated_Map = (Correlated_Map ).*(Correlated_Map > Gaussian_thresh) ;
Peaks=imregionalmax(Correlated_Map,8);
[Rows, Columns]=find(Peaks);
if max(DATA_MBonly_frame,[],'all')~=0
    SL_events_temp=repmat(struct(...
        'centroid_z',1, 'centroid_x',1,'Power',0),length(Rows),1);
    for kk = 1:length(Rows)
           centroid_loc=[-r_correc+Rows(kk)+1,-c_correc+1+Columns(kk)];
        if centroid_loc(2)<=size(DATA_MBonly_frame_high,2) && centroid_loc(1)<=size(DATA_MBonly_frame_high,1) && centroid_loc(1)>0 && centroid_loc(2)>0
            SL_events_temp(kk).centroid_z = centroid_loc(1);
            SL_events_temp(kk).centroid_x = centroid_loc(2);
            SL_events_temp(kk).Power = DATA_MBonly_frame_high(centroid_loc(1),centroid_loc(2));
        end
    end
else
    SL_events_temp=[];
end
end