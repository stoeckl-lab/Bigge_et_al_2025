%% extract contrast edges in different parts of the visual field from videofiles
% which are rotated correctly and de-shaked
% 
clear
close all

%set threshold for detecting edges
tedge=[0.01 0.025 0.05 0.1 0.5];%was 0.025 before 20240531

%% loading files
startDir='F:\Anna Backup\AG Stoeckl\Anna&Ronja\Dorsal_Ventral_Imaging\Analysis_20200930';

% Select condition and navigate to the folder
% foldername= 'Tunnel_up';
% if contains(foldername,'Tunnel_')
%     cd (['X:\ImageDifference\Analysis\TunnelAll\', num2str(foldername)]);
% else
% cd (['X:\ImageDifference\Analysis\', num2str(foldername(1:end-2)),'\', num2str(foldername)]);
% end

foldername= 'SemiCornfield_3';
if contains(foldername,'Tunnel_')
    cd ([startDir,'\TunnelAll\', num2str(foldername)]);
else
cd ([startDir,'\', num2str(foldername(1:end-2)),'\', num2str(foldername)]);
end


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
h=drawcircle('Center',c,'Radius', 575,'Drawingarea','unlimited');

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

%% extract edges from videos

%horizontal video
edges_h=extractEdge(h_video,tedge);
%vertical video
edges_v=extractEdge(v_video,tedge);

%% taking the mean

h_mean_edges=nanmean(edges_h,3);
v_mean_edges=nanmean(edges_v,3);

%% padding, to have square shape and set the real image center as center

% PADDING OF HORIZONTAL IMAGE
%desired format : 1200x1200
%actual format:
[rows,cols]=size(h_mean_edges);
% columns position:
col_diff_left=600-h_center(1);
col_diff_right=600-(cols-h_center(1));

if col_diff_left > 0 
    pad_left=nan(rows, col_diff_left);
    if col_diff_right > 0
        pad_right=nan(rows, col_diff_right);
        col_pad_edges= [pad_left,h_mean_edges, pad_right];
    else
        col_pad_edges= [pad_left,h_mean_edges(:,1:end-abs(col_diff_right))];
    end
else
    if col_diff_right > 0
        pad_right=nan(rows, col_diff_right);
         col_pad_edges= [h_mean_edges(:,abs(col_diff_left)+1:end), pad_right];
    else
        col_pad_edges= h_mean_edges(:,abs(col_diff_left)+1:end-abs(col_diff_right));
    end
end

% rows position
row_diff_up= 600-h_center(2);
row_diff_down= 600-(rows-h_center(2));

if row_diff_up > 0
    pad_up= nan(row_diff_up,1200);
    if row_diff_down > 0
        pad_down= nan(row_diff_down,1200);
        row_pad_edges= [pad_up;col_pad_edges;pad_down];
    else 
        row_pad_edges= [pad_up;col_pad_edges(1:end-abs(row_diff_down),:)];
    end
else
   if row_diff_down > 0
        pad_down= nan(row_diff_down,1200);
        row_pad_edges= [col_pad_edges(abs(row_diff_up)+1:end,:);pad_down];
   else
       row_pad_edges=col_pad_edges(abs(row_diff_up)+1:end-abs(row_diff_up),:);
   end
end

h_pad_edges= row_pad_edges;

% PADDING OF VERTICAL IMAGE
[rows,cols]=size(v_mean_edges);
% columns position:
col_diff_left=600-v_center(1);
col_diff_right=600-(cols-v_center(1));

if col_diff_left > 0 
    pad_left=nan(rows, col_diff_left);
    if col_diff_right > 0
        pad_right=nan(rows, col_diff_right);
        col_pad_edges= [pad_left,v_mean_edges, pad_right];
    else
        col_pad_edges= [pad_left,v_mean_edges(:,1:end-abs(col_diff_right))];
    end
else
    if col_diff_right > 0
        pad_right=nan(rows, col_diff_right);
         col_pad_edges= [v_mean_edges(:,col_diff_left+1:end), pad_right];
    else
        col_pad_edges= v_mean_edges(:,abs(col_diff_left)+1:end-abs(col_diff_right));
    end
end

% rows position
row_diff_up= 600-v_center(2);
row_diff_down= 600-(rows-v_center(2));

if row_diff_up > 0
    pad_up= nan(row_diff_up,1200);
    if row_diff_down > 0
        pad_down= nan(row_diff_down,1200);
        row_pad_edges= [pad_up;col_pad_edges;pad_down];
    else 
        row_pad_edges= [pad_up;col_pad_edges(1:end-abs(row_diff_down),:)];
    end
else
   if row_diff_down > 0
        pad_down= nan(row_diff_down,1200);
        row_pad_edges= [col_pad_edges(abs(row_diff_up)+1:end,:);pad_down];
   else
       row_pad_edges=col_pad_edges(abs(row_diff_up)+1:end-abs(row_diff_up),:);
   end
end

v_pad_edges=row_pad_edges;

%% merge the horizontal and vertical image
%use the maximum values 
all_edges=nanmean(cat(3,v_pad_edges,h_pad_edges),3);

%% remove the edges at the "edges" of the viewing area of hte camera
c(1)= size(all_edges,2)/2;
c(2)=size(all_edges,1)/2;
f=figure;
f.Name= 'Set the correct ROI. To continue, press any key!';
imagesc(all_edges)
h=drawcircle('Color','w','Center',c,'Radius', 560,'Drawingarea','unlimited');

disp('Set the correct ROI. Then, press any key to continue.');
pause

%use the circle object to make a mask
maskCircle=h.createMask;
maskCircle=~maskCircle;%use the invert, since we want to set all values outside of the mask to nan
close

%set all edge values outside of mask to nan
all_edges(maskCircle)=nan;

%% Calculate two hemispheres
%remember that matrix images are always flipped in y
%I keep this flip here for all
midway=size(all_edges)/2;
dorsal=false(size(all_edges));ventral=false(size(all_edges));
ventral(midway(1)+1:end,:)=true;
dorsal(1:midway(1),:)=true;

%% Calculate diagonals and quadrant coordinates for dorsal, ventral, left, right

lowerTriag=tril(all_edges);
upperTriag=triu(all_edges);

leftQ=tril(flipud(lowerTriag));
rightQ=(triu(flipud(upperTriag)));

ventralQ=flipud(triu(flipud(lowerTriag)));
dorsalQ=flipud(tril(flipud(upperTriag)));



%% Plotting all in one figure
%make a "cool" colormap
% cols=[0 0 0;68 85 196;69 137 252;42 239 160;110 254 98;165 252 69]/255;%old cols 20240531
cols=[0 0 0;68 85 196;69 137 252; 39 170 225;39 255 255;180 255 255;255 255 255]/255;
 
colmap=interp1(ceil((0:size(cols,1)-1)*(256/(size(cols,1)-1))),cols,1:256);

f= figure;
AllinOne= nan(1200,1200);
AllinOne(dorsal)= all_edges(dorsal);
AllinOne(ventral)= all_edges(ventral);

imagesc(AllinOne)
colormap(colmap)
colorbar

% hold on
% plot([0 1200],midway, 'w', 'LineWidth', 1.5);
axis equal
axis off
title(num2str(foldername))

% save figure
saveas(f,[foldername, 'contrast_heatmap.eps']);
saveas(f,[foldername, 'contrast_heatmap.fig']);
saveas(f,[foldername, 'contrast_heatmap.png']);

%% saving as mat

cdata.condition= num2str(foldername);
cdata.areaSize= sum(sum(dorsal));
cdata.indexDorsal= dorsal;
cdata.indexVentral= ventral;
cdata.MeanContrast(1)=nanmean(AllinOne(dorsal));
cdata.MeanContrast(2)=nanmean(AllinOne(ventral));
cdata.regions={'dorsal','ventral'};
cdata.MedianMag_all=AllinOne;
cdata.MeanMagnitude(1)=nanmean((dorsalQ(:)));
cdata.MeanMagnitude(2)=nanmean((ventralQ(:)));
cdata.MeanMagnitude(3)=nanmean((leftQ(:)));
cdata.MeanMagnitude(4)=nanmean((rightQ(:)));


save([foldername, 'contrast.mat'],'cdata');
disp([foldername, '.mat is saved.']);

%% extract edges
%%
function data_edges=extractEdge(v,thresh)
%EXTRACTEDGE extracts the edges in the video indicated by videoReader object v
% using the matlab inbuilt function edge
% the treshold indicates the threshold for edge detection

%set the video to first frame again
v.CurrentTime=0;

%read video data
nrFrames=round(v.Duration*v.FrameRate);

%for BW
data_edges=zeros(v.Height,v.Width,nrFrames);

%extract edges in each video frame
for u=1:nrFrames
    b = readFrame(v);
       
    b=imgaussfilt(b,4);%apply a slight blur to only extract edges that are caused by larger structures, not pixel noise
    bw=rgb2gray(b);
    %in the cornfield data the contrast is odd - problem with detecting edges
    %at the same threshold. adjust contrast therefore
    imlo=quantile(bw(:),0.01); imhi=quantile(bw(:),0.99);
    bw=imadjust(bw,double([imlo imhi])/255);

    % figure;imshow(bw);
%     [bw_edge] = edge(bw,'Prewitt',thresh);
%     figure;imshow(bw_edge);
    
%NEW 20240531
    %edge detection just returns edge, but not the strength of the
    %contrast. we could fix this by running the edge detector with multiple
    %thresholds, and then using the highest threshold for each edge, and
    %scaling entries by this threshold
    bw_edge=zeros(size(bw));
    edgeTemp_prev=zeros(size(bw));
    for t=1:length(thresh)
    edgeTemp = edge(bw,'Prewitt',thresh(end-t+1));%extract edges with highest threshold
    edgeTemp = edgeTemp-edgeTemp_prev;%remove previous edge positions from current extraction
    edgeTemp = logical(edgeTemp);
    bw_edge(edgeTemp)=thresh(end-t+1);%store threshold at the recognised edge positions
    edgeTemp_prev=edgeTemp;%store these edge positions, to then remove them from next iteration
    end

    %save the data in this matrix
    data_edges(:,:,u)=bw_edge;


end
end





