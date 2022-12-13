%% MB detection from Figure
if ~exist('skip_MB_selection')

% Select ROI
rect_roi = [];  
fig1 = figure('Position',[scr_size(1)+100 scr_size(2)+100 (scr_size(3)-200)/2 scr_size(4)-250]);
    imagesc((DATA_MBonly(:,:,1)));
    colormap(gray);
    axis square
    xlabel('Lateral (samples)')
    ylabel('Depth (samples)')
    title('Select a region of interest. For choosing all, click outside image.')
    set(gca,'FontSize',13,'Fontname','Arial')
    rect_roi = round(getrect(gca));
if  ~rect_roi(4)
    rect_roi = [1 1 size(DATA_MBonly,2)-1 size(DATA_MBonly,1)-1];  
end
close(fig1)

lateral_roi = rect_roi(3)+1;
depth_roi = rect_roi(4)+1;


% Select MBs
rect = [];  
bubble_count=0;
fig2 = figure('Position',[scr_size(1)+100 scr_size(2)+100 scr_size(3)-200 scr_size(4)-250]);
for frame = 200:no_frames;
 
    subplot(121)
    imagesc(x_axis*1e3, z_axis*1e3, (DATA_MBonly(:,:,frame)));
    colormap(gray);
    axis image
    title('Contrast-mode frame (in dB). Do not chose from this image!')
    xlabel('Lateral (mm)')
    ylabel('Depth (mm)')
    set(gca,'FontSize',13,'Fontname','Arial')

    subplot(122)
    axis square
    imagesc((DATA_MBonly(rect_roi(2)+(0:rect_roi(4)),rect_roi(1)+(0:rect_roi(3)),frame)));
    colormap(gray);
    axis square
    title('Select an isolated MB! If not, click outside image.')
    set(gca,'FontSize',13,'Fontname','Arial')
    chosen_rect = round(getrect(gca));
    if chosen_rect(4)
    bubble_count=bubble_count+1;
    rect(:,bubble_count) = [chosen_rect(1)+rect_roi(1) chosen_rect(2)+rect_roi(2) chosen_rect(3) chosen_rect(4) frame];
    end
    
    if bubble_count>= bubble_num_PSF    % if bubble_num_PSF MBs selected stop the loop
        break;
    end
end
close(fig2)
end
%% creating template

patch=zeros(patch_size,patch_size,bubble_num_PSF);
patch_r=patch;

for i=1:bubble_num_PSF
    patch(1:1+rect(4,i),1:1+rect(3,i),i)=abs(DATA_MBonly(rect(2,i)+(0:rect(4,i)),rect(1,i)+(0:rect(3,i)),rect(5,i)));
    [~,~,patch_r(:,:,i)]=icp(abs(patch(:,:,1)),patch(:,:,i));
end
patch_aver=mean(patch_r,3);

%%
x=(1:patch_size);
z=(1:patch_size);
[X,Z]=meshgrid(x,z);
tic
[Amplitude,center_x,center_z,width_x,width_z,theta]=polyGF2D(X,Z,patch_aver,noise_level_for_gaussianfit); %取对数进行二次多项式拟合
offsetp=0;
time_fit=toc;

tic
options = optimset('Display','off','TolFun',1e-4,'LargeScale','off'); %用fminunc拟合，将上面的拟合结果带入作为fminunc的初值
initpar=[Amplitude,center_x,center_z,width_x,width_z,theta,offsetp];
para=fminunc(@objfun,initpar,options,X,Z,patch_aver);
Amplitude2=para(1);center_x2=para(2);center_z2=para(3);width_x2=para(4);width_z2=para(5);theta2=para(6);offsetp2=para(7);
time_fit2=toc;

%% template for rigid registration
PSF=zeros(patch_size);
for i=1:patch_size
    for j=1:patch_size
        PSF(i,j)=(Amplitude2*exp(-(cos(0)*(j-ceil(patch_size/2))+sin(0)*(i-ceil(patch_size/2))).^2/2/width_x2^2-(-sin(0)*(j-ceil(patch_size/2))+cos(0)*(i-ceil(patch_size/2))).^2/2/width_z2^2)+0)*cos((i-(patch_size+1)/2)*2*pi/4)+...
            (Amplitude2*exp(-(cos(0)*(j-ceil(patch_size/2))+sin(0)*(i-ceil(patch_size/2))).^2/2/width_x2^2-(-sin(0)*(j-ceil(patch_size/2))+cos(0)*(i-ceil(patch_size/2))).^2/2/width_z2^2)+0)*sin((i-(patch_size+1)/2)*2*pi/4)*1i;
    end
end

figure
imagesc(abs(PSF))
colormap(gray)
save([result_path 'PSF'],'PSF')