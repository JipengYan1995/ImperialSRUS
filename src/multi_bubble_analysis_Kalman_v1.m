%%  find linked bubble ID and store in cell

tortuosity=ones(length(SL_events_linked_filtered),1)*inf;
for ii=1:length(SL_events_linked_filtered)
    if length(SL_events_linked_filtered{ii})>2
        centriod_array=[];
        for jj=1:length(SL_events_linked_filtered{ii})
            Current_centriod=[SL_events_linked_filtered{ii}(jj).centroid_x, SL_events_linked_filtered{ii}(jj).centroid_z];
            if bw_mask(Current_centriod(2),Current_centriod(1))
                centriod_array=[centriod_array;Current_centriod];
            end
        end
        
        if size(centriod_array,1)>2
           
            dis=0;
            for i=1:size(centriod_array,1)-1
                dis=dis+sqrt((centriod_array(i+1,1)-centriod_array(i,1))^2+(centriod_array(i+1,2)-centriod_array(i,2))^2);
            end
            tortuosity(ii)=dis/sqrt((centriod_array(end,1)-centriod_array(1,1))^2+(centriod_array(end,2)-centriod_array(1,2))^2);
        end
    end
    
end
tortuosity=tortuosity([tortuosity~=inf]);
tortuosity=tortuosity(~isnan(tortuosity));