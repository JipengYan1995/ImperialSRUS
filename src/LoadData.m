if strcmp(file_name(end-2:end),'dcm') || strcmp(file_name(end-2:end),'DCM')
    info = dicominfo(data_name);
    try
        frame_rate=info.CineRate;
    catch
        frame_rate = str2double(input('Input frame rate: ','s'));
    end
    video=dicomread(info,'frames',1:info.ActualFrameDuration*frame_rate);  % need to be confirmed
    
elseif strcmp(file_name(end-2:end),'avi') || strcmp(file_name(end-2:end),'AVI')
    video = VideoReader(data_name);
    try
        frame_rate=video.FrameRate;
    catch
        frame_rate = str2double(input('Input frame rate: ','s'));
    end
    
    try
        EndTime=min(EndTime,video.Duration);
    catch
        %         EndTime=min(EndTime,video.Duration);
    end
    
    try
        video=read(video,[max([round(StartTime*video.FrameRate),1]),round(EndTime*video.FrameRate)]);
    catch
        disp('Time range is out of the data \n');
    end
    %     clear
else
    disp('File format is not supported, please use avi and DICOM files')
    return
end