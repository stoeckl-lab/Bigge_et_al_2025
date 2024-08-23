%% extract optic flow in different parts of the visual field from datafiles containing
%% optic flow angles and magnitude extracted from registered videos
clear
close all
clc


%% loading files

% Select condition and navigate to the folder
% foldername= 'Tunnel_up';
% if contains(foldername,'Tunnel_')
%     cd (['X:\ImageDifference\Analysis\TunnelAll\', num2str(foldername)]);
% else
% cd (['X:\ImageDifference\Analysis\', num2str(foldername(1:end-2)),'\', num2str(foldername)]);
% end

foldername= 'Forest3_3';
if contains(foldername,'Tunnel_')
    cd (['X:\ImageDifference\Analysis_20200930\TunnelAll\', num2str(foldername)]);
else
cd (['X:\ImageDifference\Analysis_20200930\', num2str(foldername(1:end-2)),'\', num2str(foldername)]);
end


% load both AngleStacks
filename=dir('AngleStack*.mat');
if contains(filename(1).name,'turn')
    v_ang=load(filename(1).name);
else
    h_ang=load(filename(1).name);
end
if contains(filename(2).name,'turn')
   v_ang=load(filename(2).name);
else
    h_ang=load(filename(2).name);
end  

% load both magniStacks
filename=dir('MagniStack*.mat');
if contains(filename(1).name,'turn')
    v_mag=load(filename(1).name);
else
    h_mag=load(filename(1).name);

end
if contains(filename(2).name,'turn')
   v_mag=load(filename(2).name);
else
    h_mag=load(filename(2).name);
end  

h_angle=h_ang.AngleStack;
h_magni=h_mag.MagniStack;
v_angle=v_ang.AngleStack;
v_magni=v_mag.MagniStack;

% load videofiles
videofiles=dir('down_blur*.avi');
if contains(videofiles(1).name,'turn')
    v_video=VideoReader(videofiles(1).name);
else
    h_video=VideoReader(videofiles(1).name);
end
if contains(videofiles(2).name,'turn')
    v_video=VideoReader(videofiles(2).name);
else
    h_video=VideoReader(videofiles(2).name);
end

clearvars h_ang v_ang h_mag v_mag videofiles

%% Set the correct Region of Interest

h_frame=readFrame(h_video);
Img=rgb2gray(h_frame);
%c= size(Img,[2,1])/2;
c(1)= size(Img,2)/2;
c(2)=size(Img,1)/2;
f=figure;
f.Name= 'Set the correct ROI. To continue, press any key!';
imshow(Img)
h=drawcircle('Center',c,'Radius', 580,'Drawingarea','unlimited');
disp('Set the correct ROI. Then, press any key to continue.');
pause

h_center= round(h.Center);

close

v_frame=readFrame(v_video);
Img=rgb2gray(v_frame);
%c= size(Img,[2,1])/2;
c(1)= size(Img,2)/2;
c(2)=size(Img,1)/2;

f=figure;
f.Name= 'Set the correct ROI. To continue, press any key!';
imshow(Img)
h=drawcircle('Center',c,'Radius', 580,'Drawingarea','unlimited');
disp('Set the correct ROI. Then, press any key to continue.');
pause

v_center=round(h.Center);
close

%% create filtermask for both images
[hrow,hcol]=size(h_angle{1});
h_mask=createMask(hrow,hcol);

[vrow,vcol]=size(v_angle{1});
v_mask=createMask(vrow,vcol);

%% filter Angle and Magnitude with specified threshold
th= deg2rad(15);

lower_h_mask= h_mask-th;
l_ind=lower_h_mask < deg2rad(0);
lower_h_mask(l_ind)=deg2rad(360) + lower_h_mask(l_ind);
upper_h_mask= h_mask+th;
u_ind= upper_h_mask > deg2rad(360);
upper_h_mask(u_ind)= upper_h_mask (u_ind) - deg2rad(360);
% trans_ind= l_ind + u_ind;

% horizontal image
tempnan=nan(hrow,hcol);
for i=1:length(h_angle)
    h_a=h_angle{i};
    h_m=h_magni{i};
    h_ind= h_a >= lower_h_mask & h_a <= upper_h_mask;
    h_ind(l_ind)=h_a(l_ind) >= lower_h_mask(l_ind) | h_a(l_ind)<= upper_h_mask(l_ind);
    h_ind(u_ind)=h_a(u_ind) >= lower_h_mask(u_ind) | h_a(u_ind)<= upper_h_mask(u_ind);
    temp_h_a=tempnan;
    temp_h_a(h_ind)=h_a(h_ind);
    temp_h_m=zeros(hrow,hcol);
    temp_h_m(h_ind)=h_m(h_ind);
    
    h_filt_angle(:,:,i)=temp_h_a;
    h_filt_magni(:,:,i)=temp_h_m;
end

lower_v_mask= v_mask-th;
l_ind=lower_v_mask < deg2rad(0);
lower_v_mask(l_ind)=deg2rad(360) + lower_v_mask(l_ind);
upper_v_mask= v_mask+th;
u_ind= upper_v_mask > deg2rad(360);
upper_v_mask(u_ind)= upper_v_mask (u_ind) - deg2rad(360);

tempnan=nan(vrow,vcol);
for i=1:length(v_angle)
    v_a=v_angle{i};
    v_m=v_magni{i};
    v_ind= v_a >= lower_v_mask & v_a <= upper_v_mask;
    v_ind(l_ind)=v_a(l_ind) >= lower_v_mask(l_ind) | v_a(l_ind)<= upper_v_mask(l_ind);
v_ind(u_ind)=v_a(u_ind) >= lower_v_mask(u_ind) | v_a(u_ind)<= upper_v_mask(u_ind);
    temp_v_a=tempnan;
    temp_v_a(v_ind)=v_a(v_ind);
    temp_v_m=zeros(vrow,vcol);
    temp_v_m(v_ind)=v_m(v_ind);
    
    v_filt_angle(:,:,i)=temp_v_a;
    v_filt_magni(:,:,i)=temp_v_m;
end

%% taking the mean

% h_mean_angle=nanmedian(h_filt_angle,3);
% h_mean_magni=nanmedian(h_filt_magni,3);
% v_mean_angle=nanmedian(v_filt_angle,3);
% v_mean_magni=nanmedian(v_filt_magni,3);

h_mean_angle=nanmean(h_filt_angle,3);
h_mean_magni=nanmean(h_filt_magni,3);
v_mean_angle=nanmean(v_filt_angle,3);
v_mean_magni=nanmean(v_filt_magni,3);

%% padding, to have square shape and set the real image center as center

% PADDING OF HORIZONTAL IMAGE
%desired format : 1200x1200
%actual format:
h_frame=readFrame(h_video);
[rows,cols]=size(h_mean_magni);
h_frame= rgb2gray(h_frame);
% columns position:
col_diff_left=600-h_center(1);
col_diff_right=600-(cols-h_center(1));

if col_diff_left > 0 
    pad_left=nan(rows, col_diff_left);
    if col_diff_right > 0
        pad_right=nan(rows, col_diff_right);
        col_pad_magni= [pad_left,h_mean_magni, pad_right];
        col_pad_angle= [pad_left,h_mean_angle, pad_right];
        c_pad_h_frame= [pad_left, h_frame, pad_right];
    else
        col_pad_magni= [pad_left,h_mean_magni(:,1:end-abs(col_diff_right))];
        col_pad_angle= [pad_left,h_mean_angle(:,1:end-abs(col_diff_right))];
        c_pad_h_frame= [pad_left,h_frame(:,1:end-abs(col_diff_right))];
    end
else
    if col_diff_right > 0
        pad_right=nan(rows, col_diff_right);
         col_pad_magni= [h_mean_magni(:,abs(col_diff_left)+1:end), pad_right];
         col_pad_angle= [h_mean_angle(:,abs(col_diff_left)+1:end), pad_right];
         c_pad_h_frame= [h_frame(:,abs(col_diff_left)+1:end), pad_right];
    else
        col_pad_magni= h_mean_magni(:,abs(col_diff_left)+1:end-abs(col_diff_right));
        col_pad_angle= h_mean_angle(:,abs(col_diff_left)+1:end-abs(col_diff_right));
        c_pad_h_frame= h_frame(:,abs(col_diff_left)+1:end-abs(col_diff_right));
    end
end

% rows position
row_diff_up= 600-h_center(2);
row_diff_down= 600-(rows-h_center(2));

if row_diff_up > 0
    pad_up= nan(row_diff_up,1200);
    if row_diff_down > 0
        pad_down= nan(row_diff_down,1200);
        row_pad_magni= [pad_up;col_pad_magni;pad_down];
        row_pad_angle= [pad_up;col_pad_angle;pad_down];
        r_pad_h_frame= [pad_up;c_pad_h_frame;pad_down];
    else 
        row_pad_magni= [pad_up;col_pad_magni(1:end-abs(row_diff_down),:)];
        row_pad_angle= [pad_up;col_pad_angle(1:end-abs(row_diff_down),:)];
        r_pad_h_frame= [pad_up;c_pad_h_frame(1:end-abs(row_diff_down),:)];
    end
else
   if row_diff_down > 0
        pad_down= nan(row_diff_down,1200);
        row_pad_magni= [col_pad_magni(abs(row_diff_up)+1:end,:);pad_down];
        row_pad_angle= [col_pad_angle(abs(row_diff_up)+1:end,:);pad_down];
        r_pad_h_frame= [c_pad_h_frame(abs(row_diff_up)+1:end,:);pad_down];
   else
       row_pad_magni=col_pad_magni(abs(row_diff_up)+1:end-abs(row_diff_up),:);
       row_pad_angle=col_pad_angle(abs(row_diff_up)+1:end-abs(row_diff_up),:);
       r_pad_h_frame=c_pad_h_frame(abs(row_diff_up)+1:end-abs(row_diff_up),:);
   end
end

h_pad_magni= row_pad_magni;
h_pad_angle= row_pad_angle;
hori_frame= r_pad_h_frame;

% PADDING OF VERTICAL IMAGE
v_frame=readFrame(v_video);
[rows,cols]=size(v_mean_magni);
v_frame= rgb2gray(v_frame);
% columns position:
col_diff_left=600-v_center(1);
col_diff_right=600-(cols-v_center(1));

if col_diff_left > 0 
    pad_left=nan(rows, col_diff_left);
    if col_diff_right > 0
        pad_right=nan(rows, col_diff_right);
        col_pad_magni= [pad_left,v_mean_magni, pad_right];
        col_pad_angle= [pad_left,v_mean_angle, pad_right];
        c_pad_v_frame= [pad_left,v_frame, pad_right];
    else
        col_pad_magni= [pad_left,v_mean_magni(:,1:end-abs(col_diff_right))];
        col_pad_angle= [pad_left,v_mean_angle(:,1:end-abs(col_diff_right))];
        c_pad_v_frame= [pad_left,v_frame(:,1:end-abs(col_diff_right))];
    end
else
    if col_diff_right > 0
        pad_right=nan(rows, col_diff_right);
         col_pad_magni= [v_mean_magni(:,col_diff_left+1:end), pad_right];
         col_pad_angle= [v_mean_angle(:,col_diff_left+1:end), pad_right];
         c_pad_v_frame= [v_frame(:,col_diff_left+1:end), pad_right];
    else
        col_pad_magni= v_mean_magni(:,abs(col_diff_left)+1:end-abs(col_diff_right));
        col_pad_angle= v_mean_angle(:,abs(col_diff_left)+1:end-abs(col_diff_right));
        c_pad_v_frame= v_frame(:,abs(col_diff_left)+1:end-abs(col_diff_right));
    end
end

% rows position
row_diff_up= 600-v_center(2);
row_diff_down= 600-(rows-v_center(2));

if row_diff_up > 0
    pad_up= nan(row_diff_up,1200);
    if row_diff_down > 0
        pad_down= nan(row_diff_down,1200);
        row_pad_magni= [pad_up;col_pad_magni;pad_down];
        row_pad_angle= [pad_up;col_pad_angle;pad_down];
        r_pad_v_frame= [pad_up; c_pad_v_frame;pad_down];
    else 
        row_pad_magni= [pad_up;col_pad_magni(1:end-abs(row_diff_down),:)];
        row_pad_angle= [pad_up;col_pad_angle(1:end-abs(row_diff_down),:)];
        r_pad_v_frame= [pad_up; c_pad_v_frame(1:end-abs(row_diff_down),:)];
    end
else
   if row_diff_down > 0
        pad_down= nan(row_diff_down,1200);
        row_pad_magni= [col_pad_magni(abs(row_diff_up)+1:end,:);pad_down];
        row_pad_angle= [col_pad_angle(abs(row_diff_up)+1:end,:);pad_down];
        r_pad_v_frame= [ c_pad_v_frame(abs(row_diff_up)+1:end,:);pad_down];
   else
       row_pad_magni=col_pad_magni(abs(row_diff_up)+1:end-abs(row_diff_up),:);
       row_pad_angle=col_pad_angle(abs(row_diff_up)+1:end-abs(row_diff_up),:);
       r_pad_v_frame= c_pad_v_frame(abs(row_diff_up)+1:end-abs(row_diff_up),:);
   end
end

v_pad_magni=row_pad_magni;
v_pad_angle= row_pad_angle;
vert_frame= r_pad_v_frame;

%% Calculate diagonals

[~,~,left,right]=diagonalSections(h_pad_magni,'edges on');
[dorsal,ventral,~,~]=diagonalSections(v_pad_magni,'edges on');

%% Plotting all in one figure

f= figure;
AllinOne= nan(1200,1200);
AllinOne(dorsal)= v_pad_magni(dorsal);
AllinOne(ventral)= v_pad_magni(ventral);
AllinOne(left)= h_pad_magni(left);
AllinOne(right)=h_pad_magni(right);

imagesc(AllinOne)
colormap hot
colorbar

hold on
plot([0 1200],[0 1200], 'w', 'LineWidth', 1.5);
plot([0 1200],[1200 0], 'w', 'LineWidth', 1.5);
axis equal
axis off
title(num2str(foldername))

% save figure
saveas(f,[foldername, '_heatmap.eps']);
saveas(f,[foldername, '_heatmap.fig']);
saveas(f,[foldername, '_heatmap.png']);

%% plotting images
f_h= figure;
imshow(hori_frame);
f_v=figure;
imshow(vert_frame);

saveas(f_h, [foldername,'_horizontalImage.png']);
saveas(f_v,[foldername,'_verticalImage.png']);

%% saving as mat

data.condition= num2str(foldername);
data.areaSize= sum(sum(dorsal));
data.indexDorsal= dorsal;
data.indexVentral= ventral;
data.indexLeft=left;
data.indexRight=right;
data.MeanMagnitude(1)=nanmean(AllinOne(dorsal));
data.MeanMagnitude(2)=nanmean(AllinOne(ventral));
data.MeanMagnitude(3)=nanmean(AllinOne(left));
data.MeanMagnitude(4)=nanmean(AllinOne(right));
data.regions={'dorsal','ventral','left','right'};
data.MedianMag_all=AllinOne;

save([foldername, '.mat'],'data');
disp([foldername, '.mat is saved.']);

%% create Mask function
%%
function mask=createMask(row,col)
%CREATEMASK makes the angle mask to filter only translational optic flow
% it uses the size of the input matrix to create a new mask with correct
% angles in degree

center= [row,col]/2;
zero_deg= [center(1),col];

vec_center=[center(1)-zero_deg(1), center(2)- zero_deg(2)];
mask=zeros(row,col);

for r=1:row
    for c=1:col
        p= [r,c];
        vec_p= [center(1)-p(1), center(2)-p(2)];
        alpha= acosd(dot(vec_center,vec_p)/(norm(vec_center)*norm(vec_p)));
        beta= 360-alpha;
        
        if r<=center(1)
            angle=beta;
        else
            angle=alpha;
        end
        
        mask(r,c)=deg2rad(angle);
    end
end
end





