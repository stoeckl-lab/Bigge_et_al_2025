function autoTrack_tunnel_dorsal(thresh1,thresh2, windowSize)
%load videos of the tunnel and digitise the moths' positions


lengthTunnel=910;%length of tunnel in mm (visible portion of the tunnel)
widthTunnel=300;%length of tunnel in mm

%check if threshold parameters were set, if not, use default
if nargin<2
    thresh1=20;%threshold for deleting image differences before Gaussian filtering
    thresh2=1;%second threshold  for deleting image differences after Gaussian filtering (needs to be lower because smoothing reduces peak size of noise and signal)
    windowSize=60;%window around the moth - everything outside of this window is deleted so that automatic tracking works better. adjust with size of camera and speed of animal
elseif nargin==2 && nargin <3
    windowSize=60;
end

%% read existing tunnel file
if exist('tunnel_pixels.mat')~=0
    if input('load existing tunnel file? Press 1')==1
        load('tunnel_pixels.mat')
        
        preRotTunnelPxs=tunnelPxs;
        %determine how the tunnel is rotated with respect to the x axis
        %use upper and lower bound of tunnel and make average
        alpha1=atand((tunnelPxs(1,2)-tunnelPxs(2,2))/(tunnelPxs(1,1)-tunnelPxs(2,1)));
        alpha2=atand((tunnelPxs(3,2)-tunnelPxs(4,2))/(tunnelPxs(3,1)-tunnelPxs(4,1)));
        alpha=(alpha1+alpha2)/2;
        
        %rotation matrix
        R=[cosd(-alpha) -sind(-alpha); sind(-alpha) cosd(-alpha)];
        %rotate all coordinates
        tunnelPxs=[R*tunnelPxs']';
    end
end

%% read file names
%to reanalyse data that has already been tracked (and where we know the
%tracker works, can just select the filenames of these data files, and
%retrack them without checking all the images
% filenames=uigetfile('.mat','multiselect','on');
% for i=1:length(filenames)
%     filenames{i}(end-3:end)='.mp4';
%     filenames{i}(1:10)=[];
% end
% if iscell(filenames)==0
%     temp=filenames; clear filenames;
%     filenames=cell(1,1);
%     filenames{1}=temp; clear temp;
% end

% %this is to select filenames in the "normal" mode from videos
filenames=uigetfile('.mp4','multiselect','on');
if iscell(filenames)==0
    temp=filenames; clear filenames;
    filenames=cell(1,1);
    filenames{1}=temp; clear temp;
end

%% analysis of each video
for j=1:length(filenames)
    disp(filenames{j});
    
    v = VideoReader(filenames{j});
    
    nrFrames=round(v.Duration*v.FrameRate);
    dataBW_all=zeros(v.Height,v.Width,2);%just save the current and previous frame
    dataBW_diff=zeros(v.Height,v.Width,nrFrames-1);
    
    coordinates=nan(nrFrames,4);
    areas=nan(nrFrames,3);
    
    count=0;
    while hasFrame(v)
        count=count+1;
        raw = readFrame(v);
        dataBW=raw(:,:,1);%to get rid of red patterns
%         dataBW=rgb2gray(raw);
        
        %flip the yaxis (rows) of the image, because when extracting coordinates from
        %pictures into matrices, yaxis becomes flipped
        %therefore, flip pic to have the matrix later correct compared to the original
        %imshow and imagesc corrects this, but not corrected in coordinate
        %system of the tracking coordinates
        dataBW=dataBW(end:-1:1,:);
        if count ==1
            exampleFrame=dataBW;
        end
        
        %get the tunnel dimensions
        if j==1 && count ==1 && exist('tunnelPxs')==0
            figure;imshow(dataBW,[0 220]);
            disp('mark the tunnel dimensions, upper border: left - right, lower border: left - right');
            [tunnelPxs]=ginput(4);
            save('tunnel_pixels.mat','tunnelPxs');
            preRotTunnelPxs=tunnelPxs;
            
            %determine how the tunnel is rotated with respect to the x axis
            %use upper and lower bound of tunnel and make average
            alpha1=atand((tunnelPxs(1,2)-tunnelPxs(2,2))/(tunnelPxs(1,1)-tunnelPxs(2,1)));
            alpha2=atand((tunnelPxs(3,2)-tunnelPxs(4,2))/(tunnelPxs(3,1)-tunnelPxs(4,1)));
            alpha=(alpha1+alpha2)/2;
            
            %rotation matrix
            R=[cosd(-alpha) -sind(-alpha); sind(-alpha) cosd(-alpha)];
            %rotate all coordinates
            tunnelPxs=[R*tunnelPxs']';
        end
        
        %rotate image - checked it rotates correctly
        rotImage=imrotate(dataBW,alpha,'bicubic','crop');
        
        %set everything outside the tunnel dimensions to 0
        rotImage(1:round(nanmin(tunnelPxs(1:2,2)))-100,:)=0;
        rotImage(round(nanmax(tunnelPxs(3:4,2)))+100:end,:)=0;
        
        if count==1
            dataBW_all(:,:,1)=rotImage; %save the first image for background subtraction
            
            %generate image differences
        elseif count>1
            
            dataBW_all(:,:,2)=rotImage; %store the image that is compared to the start image here
            
            %do differences against first frame, but needs
            %to be one that is "moth free"
            temp_diff=abs(dataBW_all(:,:,2)-dataBW_all(:,:,1));
            %image differences when animal present are quite big. discard smaller
            %changes due to fluctuations on the sensor
            temp_diff(temp_diff<thresh1)=0;
%             dataBW_diff(:,:,count-1)=temp_diff;
            
            %here is the tracking
            labelingImage=imgaussfilt(temp_diff,1)>thresh2;%gaussian filter image, to get rid of noise
%             labelingImage=temp_diff>thresh2;%gaussian filter image, to get rid of noise

            %if there is a moth detected, set everything outside a square around the detected moth to 0
            %the square needs to be big enough to allow the detected moth to move to this from
            %the previous tracking point
            %even if there wasn't a detection for one frame, use the
            %coordinate from the previous one to assume that is where the
            %moth is
            if count>2 && sum(isnan(coordinates(1:count-1,1))~=1)~=0
                lastCoordInd=find(isnan(coordinates(1:count-1,1))==0,1,'last');
                labelingImage(1:round(coordinates(lastCoordInd,2))-windowSize,:)=0;
                labelingImage(round(coordinates(lastCoordInd,2))+windowSize:end,:)=0;
                labelingImage(:,1:round(coordinates(lastCoordInd,1))-windowSize)=0;
                labelingImage(:,round(coordinates(lastCoordInd,1))+windowSize:end)=0;
            end
            
            dataBW_diff(:,:,count-1)=labelingImage;

            [labeledImage, numSpots] = bwlabel(labelingImage);
            %this allows for multiple spots, and picks the largest one, so
            %that little changes in background do not affect the tracking
            %it does not work to exclude other moths in the tunnel
            if numSpots>=1 
                props = regionprops(labeledImage, 'Centroid', 'Area','MajorAxisLength','MinorAxisLength');
                if nanmax([props.Area])>25 %make sure not to track too small points
                if numel([props.Area])>1
                    if sum(isnan(coordinates(1:count-1,1))==0)==0 %so if this is the first entry,
                        %                         then the marker has to be at the edges of the tunnel
                        tempCoord=reshape([props.Centroid],numel([props.Centroid])/2,2)';
                        entryWindow=max(tunnelPxs(:,1))*(windowSize/1000);%make start and end the same window size
                        props(tempCoord(:,1)>entryWindow & tempCoord(:,1)<max(tunnelPxs(:,1))-entryWindow)=[];
                        %then if there are still more than one, take the
                        %largest
                        if ~isempty(props)
                         if numel([props(max([props.Area])==[props.Area]).Area])==1
                        coordinates(count,1:2)=props(max([props.Area])==[props.Area]).Centroid;
                        areas(count,:)=[props(max([props.Area])==[props.Area]).Area props(max([props.Area])==[props.Area]).MajorAxisLength props(max([props.Area])==[props.Area]).MinorAxisLength];
                        else 
                            coordinates(count,1:2)=[nan nan];
                            areas(count,:)=[nan nan nan];
                         end
                        end
                    else %take the largest area
                        if numel([props(max([props.Area])==[props.Area]).Area])==1
                        coordinates(count,1:2)=props(max([props.Area])==[props.Area]).Centroid;
                        areas(count,:)=[props(max([props.Area])==[props.Area]).Area props(max([props.Area])==[props.Area]).MajorAxisLength props(max([props.Area])==[props.Area]).MinorAxisLength];
                        else 
                            coordinates(count,1:2)=[nan nan];
                            areas(count,:)=[nan nan nan];
                        end
                    end
                elseif numel([props.Area])==1
                    %so if this is the first entry, then the marker has to be at the edges of the tunnel
                    if sum(isnan(coordinates(1:count-1,1))==0)==0
                        tempCoord=reshape([props.Centroid],numel([props.Centroid])/2,2);
                        entryWindow=max(tunnelPxs(:,1))*(windowSize/1000);%make start and end the same window size
                        props(tempCoord(:,1)>entryWindow & tempCoord(:,1)<max(tunnelPxs(:,1))-entryWindow)=[];
                        %save moth coordinates
                        if isempty(props)==0
                            coordinates(count,1:2)=props.Centroid;
                            areas(count,:)=[props.Area props.MajorAxisLength props.MinorAxisLength];
                        end
                    else
                        if isempty(props)==0
                            coordinates(count,1:2)=props.Centroid;
                            areas(count,:)=[props.Area props.MajorAxisLength props.MinorAxisLength];
                        end
                    end
                end
%                 figure;imagesc(labelingImage);hold on;plot(coordinates(count,1),coordinates(count,2),'o');
                end
            end
        end
    end
    
    
    %% save tunnel coordinates
    coordinates(1:4,3:4)=tunnelPxs;
    
    % assume field of view is 1m (length of small tunnel)
    % can either use tunnel dimensions at set length of tunnel (currently
    % points indicate left and right extremes)
    % or just use pixel dimensions of camera...
    ymmperpx=widthTunnel/(diff(tunnelPxs([1,3],2))+diff(tunnelPxs([2,4],2)))*2;
    xmmperpx=lengthTunnel/nanmean(abs(diff(tunnelPxs(:,1))));
    mmperpx=(xmmperpx+ymmperpx)/2;
    
    moth=mmperpx*coordinates(:,1:2);
    tunnel=mmperpx*tunnelPxs;
    
    
    %% plot result and save each coordinates
    exampleFrame=imrotate(exampleFrame,alpha,'bicubic','crop');
    tunnelCenter=[[tunnelPxs(1,1), tunnelPxs(1,2)-0.5*(tunnelPxs(1,2) - tunnelPxs(3,2))];[tunnelPxs(2,1), tunnelPxs(1,2)-0.5*(tunnelPxs(1,2) - tunnelPxs(3,2))]];
    
    
    %plot real tunnel with tracks overlaid
    f1=figure('Position',[263 484 1600 420]);
    hold on;
    imshow(exampleFrame);hold on;
    green = cat(3, zeros(size(exampleFrame)),ones(size(exampleFrame)), ones(size(exampleFrame)));
    h = imshow(green);
    set(h, 'AlphaData', nansum(dataBW_diff,3))
    plot(tunnelPxs(1:2,1),tunnelPxs(1:2,2),'w*-');
    plot(tunnelPxs(3:4,1),tunnelPxs(3:4,2),'w*-');
    plot(tunnelCenter(:,1),tunnelCenter(:,2),'w*--');
    plot(coordinates(:,1),coordinates(:,2),'r-*');
    
    xlim([floor(tunnelPxs(1,1))-2 ceil(tunnelPxs(2,1))+2])
    ylim([floor(tunnelPxs(1,2))-20 ceil(tunnelPxs(3,2))+20])
    set(gca,'ydir','normal');
    set(gca,'PlotBoxAspectRatio',[3.33 1 1])
    title(filenames{j}(1:end-4), 'Interpreter', 'none');
    box on
    
    f2=figure('Position',[263 484 1600 420]);
    imagesc(nansum(dataBW_diff,3),[min(dataBW_diff(:)), 2.5*max(dataBW_diff(:))]);hold on;
    
    plot(tunnelPxs(1:2,1),tunnelPxs(1:2,2),'w*-');
    plot(tunnelPxs(3:4,1),tunnelPxs(3:4,2),'w*-');
    plot(tunnelCenter(:,1),tunnelCenter(:,2),'w*--');
    plot(coordinates(:,1),coordinates(:,2),'r-*');
    
    xlim([floor(tunnelPxs(1,1))-2 ceil(tunnelPxs(2,1))+2])
    ylim([floor(tunnelPxs(1,2))-20 ceil(tunnelPxs(3,2))+20])
    set(gca,'ydir','normal');
    set(gca,'PlotBoxAspectRatio',[3.33 1 1])
    title(filenames{j}(1:end-4), 'Interpreter', 'none');
    box on
    %
    pause(1) %somehow need a short break before getting input, or the image might not be drawn before input is requested,
    %and then the script gets stuck

    evaluation=input('press 1 for saving, press 2 for manual tracking, 3 for removing intruder moth and 0 for not saving');
    if evaluation==1
        save(['auto_data_',filenames{j}(1:end-4)],'coordinates','moth','tunnel','preRotTunnelPxs','areas')
    elseif evaluation==2
        DLTdv6_tunnel
        pause
        
    elseif evaluation==3
        %% try to delete a "second moth" in the tunnel, and see if that makes it
        %better
        % delete coordinates that are not in direct contact with the longest
        % continuous trace by counting the number of nan~= entries. allow for
        % some mistakes in tracking, so that 2 nan entries behind each other
        % are fine
        nanInds=isnan(coordinates(:,1))==0;
        startInd=find(diff(nanInds)==1);endInd=find(diff(nanInds)==-1);
        if numel(endInd)<numel(startInd)
            endInd(end+1)=length(nanInds);
        elseif numel(endInd)>numel(startInd)
            startInd=[1; startInd];
        end
        lengthTrack=zeros(length(startInd),1);
        for l=1:length(startInd)
            lengthTrack(l)=abs(coordinates(startInd(l)+1,1)-coordinates(endInd(l),1)-1);
        end
        %turn indices for longest track 0, keep the others 1 so they can be
        %used for deleting
        nanInds(startInd(lengthTrack==max(lengthTrack)):endInd(lengthTrack==max(lengthTrack)))=0;
        coordinates(nanInds,1:2)=nan;
        moth=mmperpx*coordinates(:,1:2);%resave coordinates for moth
        areas(nanInds,1:3)=nan;
        
        %plot updated coordinates again
        f2=figure('Position',[263 484 1600 420]);
        imagesc(nansum(dataBW_diff,3),[min(dataBW_diff(:)), 2.5*max(dataBW_diff(:))]);hold on;
        
        plot(tunnelPxs(1:2,1),tunnelPxs(1:2,2),'w*-');
        plot(tunnelPxs(3:4,1),tunnelPxs(3:4,2),'w*-');
        plot(tunnelCenter(:,1),tunnelCenter(:,2),'w*--');
        plot(coordinates(:,1),coordinates(:,2),'r-*');
        
        xlim([floor(tunnelPxs(1,1))-2 ceil(tunnelPxs(2,1))+2])
        ylim([floor(tunnelPxs(1,2))-20 ceil(tunnelPxs(3,2))+20])
        set(gca,'ydir','normal');
        set(gca,'PlotBoxAspectRatio',[3.33 1 1])
        title(filenames{j}(1:end-4), 'Interpreter', 'none');
        box on
        %
        %
        pause(1)
        evaluation=input('did this fix the problem? press 1 for saving, 2 for manual tracking, and 0 for not saving');
        if evaluation==1
            save(['auto_data_',filenames{j}(1:end-4)],'coordinates','moth','tunnel','preRotTunnelPxs','areas')
        elseif evaluation==2
            DLTdv6_tunnel
            pause
        end
        
    end
    %
    close(f1,f2)
    
end

end
