function analyse_tunnel_tracks(varargin)
%quantitative analysis of individual hawkmoth tracks through the tunnel

%tunnel limits are an essential analysis parameter. allow user input, so
%that it can be consistent throughout one project
if isempty(find(strcmp(varargin,'tunnelLimits')))==0
    tunnelLimits = varargin{find(strcmp(varargin,'tunnelLimits'))+1};
    fprintf('tunnel limits for analysis: %u - %u',tunnelLimits(1),tunnelLimits(2));
else
    tunnelLimits=[200 800];%for dorsal ventral analysis
    fprintf('tunnel limits for analysis: %u - %u',tunnelLimits(1),tunnelLimits(2));
end


%load filenams
files=uigetfile('multiselect','on','*data*.mat');
if iscell(files)==0
    temp=files;
    clear files
    files=cell(1,1); files{1}=temp;
end

Fs=50;%sampling rate Hz




%initialise
n=length(files);
% perc_half=zeros(n,2);
all_speed=cell(n,1);
all_pos=cell(n,1);
var_pos=nan(n,1);
speeds=nan(n,3);
areas=nan(n,1);
all_angles=cell(n,1);
avg_angles=nan(n,1);
tortuosity=nan(n,1);
positions=nan(n,5);
speedBins=0:100:2500;
histSpeed=nan(n,length(speedBins)-1);
allTracks=cell(n,1);
ID=nan(n,1);%identity of each animal. some experiments have this in the name of the file

bin_tunnel=[0:100:1000]; %bin tunnel into 5cm bins
bin_speed=nan(n,length(bin_tunnel)-1);

for i=1:length(files)
    load(files{i});
    if strfind(files{i},'data')~=0 && ~isempty(moth) %%added 20240213
        load(files{i});
        %         disp(files{i});
        %sort into left and right, so that all start from the same side
        %(the left side)
        if strfind(files{i},'data_R')
            moth(:,1)=1000- moth(:,1);
        end

        %in the beginnging, autoTrack tunnel did not properly delete
        %intruders from the "moth" and "areas" variable, just from coordinates - thus
        %use coordinates to clean up the "moth" variable
        if size(moth,1)==size(coordinates,1)
            moth(isnan(coordinates(:,1)),:)=nan;
            areas(isnan(coordinates(:,1)),:)=nan;
        end

        %sometimes tracks are from the right, even though it says left moth
        %- because multiple tracks in one video, and the one that got
        %selected was from the left
        %add the last part because sometimes a moth comes in after the first one
        %that starts "more left" and then looks like the moth is going
        %right, when it doesn't
        if moth(find(~isnan(moth(:,1)),1,'first'),1) > moth(find(~isnan(moth(:,1)),1,'last'),1) && moth(find(~isnan(moth(:,1)),1,'first'),1)>500
            moth(:,1)=1000- moth(:,1);
        end

        %remove points from areas that are not used for tunnel
        %analysis - generate indices, also for when moth is out of
        %limit and then flies back in (mainly at exit)
        startInd=find(moth(:,1)<tunnelLimits(1),1,'last');
        endInd=find(moth(:,1)>tunnelLimits(2),1,'first');

        if ~isempty(startInd) && ~isempty(endInd)
            %if an intruder moth came in later, it might confuse this the start
            %and end point
            if startInd>endInd
                startInd=find(moth(1:endInd,1)<tunnelLimits(1),1,'last');
            end

            %get index of individual moth
            ind1=strfind(files{i},'L');ind2=strfind(files{i},'R');
            if isempty(ind1)==0
                uscoreInd=strfind(files{i}(ind1:end),'_');
                if uscoreInd(1)>2
                    if ~isempty(str2num([files{i}(ind1+1:uscoreInd(1)-2+ind1)]))
                        ID(i)=str2num([files{i}(ind1+1:uscoreInd(1)-2+ind1)]);
                    end
                end
            elseif isempty(ind2)==0
                uscoreInd=strfind(files{i}(ind2:end),'_');
                if uscoreInd(1)>2
                    if ~isempty(str2num([files{i}(ind2+1:uscoreInd(1)-2+ind2)])) %to make sure it doesn't convert letter identifiers for tunnel conditions
                        ID(i)=str2num([files{i}(ind2+1:uscoreInd(1)-2+ind2)]);
                    end
                end
            end


            %analyse entire tunnel for tunnel wide statistics
            [incr_dist_temp, dummy]=euclid_dist(moth);
            for j=1:length(bin_tunnel)-1


                %calculate speed
                if sum(isnan(find(moth(1:end-1,1)<bin_tunnel(j+1)& moth(1:end-1,1)>bin_tunnel(j))))
                else
                    bin_speed(i,j)=nanmean(incr_dist_temp(find(moth(1:end-1,1)<bin_tunnel(j+1)& moth(1:end-1,1)>bin_tunnel(j)))*Fs);
                end
            end

            %some traces have the upper tunnel border up, and the lower one down.
            %have to swap them! this is the ones that were analysed with
            %DLTdv6!!!
            if tunnel(3,2)<tunnel(1,2)
                temp=tunnel(1:4,2);tunnel(3:4,2)=temp(1:2);tunnel(1:2,2)=temp(3:4);
            end

            tunnelCenter=[[tunnel(1,1), tunnel(1,2)-0.5*(tunnel(1,2) - tunnel(3,2))];[tunnel(2,1), tunnel(1,2)-0.5*(tunnel(1,2) - tunnel(3,2))]];


            %select only forward flights, no turning in middle or and flying back (even
            %as part of flight)
            %         forwards=diff(moth(:,1))<0;
            %         if sum(forwards)<60

            %collect all tracks for paper source files
            flightTrack=[moth(:,1) moth(:,2)-nanmean(tunnelCenter(:,2))];
            flightTrack(isnan(flightTrack(:,1)),:)=[];
            allTracks{i}=flightTrack;

            %extract area of tracking point for size estimation
            if exist('areas') && nansum(nansum(~isnan(areas).*(areas~=0)))>0
                %remove points from areas that are not used for tunnel
                %analysis
                temp_area=areas;%indDel=(moth(:,1)<tunnelLimits(1))+(moth(:,1)>tunnelLimits(2))>0;
                temp_area(endInd:end,:)=[];temp_area(1:startInd,:)=[];
            else
                all_areas{i}=nan;
            end

            %analyse only center part of the tunnel
            moth(endInd:end,:)=[];moth(1:startInd,:)=[];

            %make sure that areas has same size as moth
            if exist('areas') && nansum(nansum(~isnan(areas).*(areas~=0)))>0
                %remove points from areas that are not used for tunnel
                %analysis
                temp_area(isnan(moth(:,1)),:)=[];
                all_areas{i}=temp_area;
            end

            moth(isnan(moth(:,1)),:)=[];

            if ~isempty(moth)

                %calculate flight speed
                [incr_dist, all_angles{i}]=euclid_dist(moth);
                all_speed{i}=incr_dist*Fs;
                speeds(i,1)=nanmedian(all_speed{i});

                speeds(i,2)=quantile(all_speed{i},.90);%compute the max. 10% speed
                speeds(i,3)=nanstd(all_speed{i});%use standard deviation because more robust to sample size

                %generate a histogram of speeds, to show portions of hovering and
                %straight flight
                [histSpeed(i,:),histSpeed_edges] = histcounts(all_speed{i},speedBins);
                %make figure of flight track and speed histogram
                %             figure('Position',[400 400 1200 300]);
                %             subplot(1,3,[1,2]); plot(moth(:,1),moth(:,2),'o-');
                %             title(files{i}, 'Interpreter', 'none')
                %             ylim([tunnel(1,2) tunnel(3,2)]);
                %
                %             subplot(1,3,3),bar(histSpeed_edges(1:end-1)+0.5*nanmean(diff(histSpeed_edges)),histSpeed(i,:))
                %             ylim([0,15]);xlim([0 2500]);
                %             title(sprintf('avg speed: %3.0f',speeds(i,1)))

                % calculate angles between consecutive tracking points
                avg_angles(i)=abs(nanmean(all_angles{i}));

                %tortuosity
                %smooth tracks a bit to take out effect of little wobbles and ripples, and especially hovering on spot
                smoothMoth=moth;smoothMoth(:,2)=smooth(moth(:,2),5);
                [tort_dist, ~]=euclid_dist(smoothMoth);
                %could calculate tortuosity for real "bird distance" or path, or for
                %distance of tunnel, considering the moth should traverse it in a straight
                %line from one side to the other ("bird distance":
                %euclid_dist(moth(1,:),moth(end,:)
                
                tortuosity(i)=nansum(tort_dist)/euclid_dist(smoothMoth(1,:),[smoothMoth(end,1) smoothMoth(end,2)]);

                %position of the moths in the tunnel
                all_pos{i}=moth(:,2)-nanmean(tunnelCenter(:,2));
                [histPosition(i,:),histPos_edges] = histcounts(all_pos{i},[-150:25:150],'Normalization','probability');

                positions(i,1)=nanmedian(all_pos{i});%median position !! This weighs by how much time the moth spend where!!
                %but actually, if we want to use median position as a
                %descriptor of the path shape and thus where the animal flew,
                %we need to interpolate the path, and then take the median
                %position! because otherwise when it flew fast these "points"
                %are underrepresented
                if length(moth(:,1))==length(unique(moth(:,1))) && length(moth(:,2))==length(unique(moth(:,2)))
                    intTrace=interp1(moth(:,1),moth(:,2)-nanmean(tunnelCenter(:,2)),[tunnelLimits(1):20:tunnelLimits(2)],[],'extrap');%interpolate to that outside range they are nan
                else
                    intTrace=nan;
                end
                positions(i,4)=nanmedian(intTrace);%variation (s.d.) in position
                %but 2D interpolation currently isn't working  for animals moving backwards, so I made the step super coarse!!
                %still some glitches for very loopy trajectories but better on
                %the whole

                positions(i,2)=nanstd(all_pos{i});%variation (s.d.) in position

                %crossover index (did moth cross from one half to the other
                %half)
                first=nanmean(all_pos{i}(moth(:,1)<300));second=nanmean(all_pos{i}(moth(:,1)>700));
                positions(i,3)=first-second;

                %position in first half of the tunnel. in case the moth just
                %quickly reacted to the pattern, but then flew back to the
                %middle
                positions(i,6)=nanmedian(moth(1:floor(size(moth,1)/2),2))-nanmean(tunnelCenter(:,2));

                %crossover position in first half of the tunnel. in case the moth just
                %quickly reacted to the pattern, but then flew back to the
                %middle
                first=nanmean(all_pos{i}(moth(:,1)<350));second=nanmean(all_pos{i}(moth(:,1)>350));
                positions(i,7)=first-second;

                %percentage lateral travel
                xtravel=nansum(abs(diff(moth(:,1))));
                ytravel=nansum(abs(diff(moth(:,2))));
                positions(i,5)=ytravel/xtravel;
                %alternative frame by frame analysis 20240220
%                             positions(i,5)=nanmean(abs(diff(moth(:,2))./diff(moth(:,1))));

                if exist('area')==1
                    area(area==0)=[];
                    areas(i)=nanmean(area);
                    clear area
                end

                %check6
                % figure;plot(smoothMoth(:,1),smoothMoth(:,2),'b.-');hold on;plot(smoothMoth([1,end],1),[smoothMoth(1,2);smoothMoth(1,2)],'r.-');ylim([tunnel(1,2), tunnel(4,2)])
                % plot(tunnelCenter(:,1),tunnelCenter(:,2),'k--');xlim(tunnelLimits);
                % disp(speeds(i,:))
                %         else
                %         end
            end
        end
    end
end



%% PLOTS

%plot speed along tunnel
x=bin_tunnel(2:end)-nanmean(diff(bin_tunnel))/2;
avg=nanmean(bin_speed);
err=nanstd(bin_speed)/(sqrt(size(bin_speed,1))-1);
% figure;hold on;
% plot(x,bin_speed','color',[0.5 0.5 0.5])
% errorbar(x,avg,err,'ko-','MarkerFaceColor','k')
% xlabel('pos in tunnel [mm]');ylabel('speed [mm/s]');

%summary stats
figure('position',[400   300   1000   600]); hold on;
xmin=0.8;xmax=1.2;spread=0.02;
dotCol=[0.5 0.5 0.5];
boxCol=[0 0 0];

subplot(2,4,1);hold on;
plot(spread*randn(n,1)+ones(n,1),positions(:,1),'.k','MarkerSize',16,'color',dotCol);
% plot(0.05*randn(n,1)+ones(n,1),perc_half(:,2),'.r','MarkerSize',16);
b=boxplot(positions(:,1),'color',boxCol);set(b,{'linew'},{2})
plot([0.5 1.5],[0 0 ],'k--');
title('median pos');ylim([-150 150]);xlim([xmin xmax]);
ylabel('mm');

subplot(2,4,2);hold on;
plot(spread*randn(n,1)+ones(n,1),speeds(:,1),'.k','MarkerSize',16,'color',dotCol);
b=boxplot(speeds(:,1),'color',boxCol);set(b,{'linew'},{2})
title('avg speed');ylim([0 2500]);xlim([xmin xmax]);
ylabel('mm/s');

tunnelCenter=[[tunnel(1,1), tunnel(1,2)-0.5*(tunnel(1,2) - tunnel(3,2))];[tunnel(2,1), tunnel(1,2)-0.5*(tunnel(1,2) - tunnel(3,2))]];

subplot(2,4,3);hold on;
plot(spread*randn(n,1)+ones(n,1),positions(:,2),'.k','MarkerSize',16,'color',dotCol);
% plot(0.05*randn(n,1)+ones(n,1),perc_half(:,2),'.r','MarkerSize',16);
b=boxplot(positions(:,2),'color',boxCol);set(b,{'linew'},{2})
title('pos variation');ylim([0 100]);xlim([xmin xmax]);
ylabel('std (mm)');

subplot(2,4,4);hold on;
plot(spread*randn(n,1)+ones(n,1),tortuosity,'.k','MarkerSize',16,'color',dotCol);
b=boxplot(tortuosity,'color',boxCol);set(b,{'linew'},{2})
title('avg tortuosity');ylim([1 1.5]);xlim([xmin xmax]);
ylabel('length/distance');

subplot(2,4,5);hold on;
plot(spread*randn(n,1)+ones(n,1),positions(:,3),'.k','MarkerSize',16,'color',dotCol);
b=boxplot(positions(:,3),'color',boxCol);set(b,{'linew'},{2})
title('crossover index');ylim([-150 150]);xlim([xmin xmax]);
ylabel('mm');

subplot(2,4,6);hold on;
plot(spread*randn(n,1)+ones(n,1),positions(:,5),'.k','MarkerSize',16,'color',dotCol);
% plot(0.05*randn(n,1)+ones(n,1),perc_half(:,2),'.r','MarkerSize',16);
b=boxplot(positions(:,5),'color',boxCol);set(b,{'linew'},{2})
title('prop lateral mov');ylim([0 2]);xlim([xmin xmax]);
ylabel('y/x travel');

subplot(2,4,7);hold on;
plot(spread*randn(n,1)+ones(n,1),positions(:,6),'.k','MarkerSize',16,'color',dotCol);
b=boxplot(positions(:,6),'color',boxCol);set(b,{'linew'},{2})
plot([0.5 1.5],[0 0 ],'k--');
title('pos 1st half');ylim([-150 150]);xlim([xmin xmax]);
ylabel('mm');

subplot(2,4,8);hold on
plot(spread*randn(n,1)+ones(n,1),positions(:,7),'.k','MarkerSize',16,'color',dotCol);
b=boxplot(positions(:,7),'color',boxCol);set(b,{'linew'},{2})
title('crossover index 1stH');ylim([-150 150]);xlim([xmin xmax]);
ylabel('mm');



%% analyse speed
histSpeedS=struct;
histSpeedS.counts=histSpeed;
histSpeedS.edges=histSpeed_edges;

% plot position histogram
[histPosition(i,:),histPos_edges] = histcounts(all_pos{i},[-150:25:150]);

% figure;shadedErrorBar_anna(histSpeed_edges(1:end-1)+0.5*nanmean(diff(histSpeed_edges)),nanmean(histSpeed,1),2*nanstd(histSpeed,1)/(sqrt(size(histSpeed,1))));
% figure;hold on;
% bar(histPos_edges(1:end-1)+0.5*nanmean(diff(histPos_edges)),nanmean(histPosition,1));
% errorbar(histPos_edges(1:end-1)+0.5*nanmean(diff(histPos_edges)),nanmean(histPosition,1),2*nanstd(histPosition,1)/(sqrt(size(histPosition,1))),'k.');
% xlim([-150 150]);

% test correlation speed and tunnel position
% pos=abs(positions(:,1));%make this absolute for symmetrical tunnel patterns, because there side shouldn't matter
% speed=speeds(:,1);%1 is avg speed, 2 is max speed
% pos(isnan(pos))=[];
% speed(isnan(speed))=[];
% 
% [R,P] = corrcoef(pos,speed,'rows','pairwise');
% 
% p = polyfit(pos,speed,1);
% f = polyval(p,pos);
% 
% pos(isnan(pos))=[];speed(isnan(pos))=[];
% pos(isnan(speed))=[];speed(isnan(speed))=[];
% 
% figure;plot(pos,speed,'ob');ylabel('speed');xlabel('position')
% hold on
% plot(pos,f,'k')
% xlim([0 150])
% title(sprintf('R=%1.1f, p=%1.3f',R(2),P(2)));

%% save
c=cd;
%get rid of underscore as they show up nasitly
% indScore=strfind(c,'_');c(indScore)=' ';
% indc=strfind(c,'2018');
% suptitle(c(indc:end));

%read out condition
currPath=pwd;inds=strfind(currPath,'\');
condition=currPath(inds(end)+1:end);

save('comp_measures.mat','all_speed','tortuosity','all_pos','speeds','positions','histSpeedS','ID','condition','all_areas','allTracks')

%save for publication source data
% position=positions(isnan(positions(:,1))==0,1);inflightvar=positions(isnan(positions(:,1))==0,2);speed=speeds(isnan(positions(:,1))==0,1);max_speed=speeds(isnan(positions(:,1))==0,2);
allTracks(isnan(positions(:,1)))=[];
save('all_tracks.mat','allTracks');
% save('summary_data.mat','position','inflightvar','speed','max_speed');

%% plot individual points with colours for each individual animal
if sum(isnan(ID))==0
    %summary stats
    figure('position',[400   300   1000   600]); hold on;
    xmin=0.8;xmax=1.2;spread=0.02;
    dotCol=[0.5 0.5 0.5];
    boxCol=[0 0 0];
    col=hsv;
    individuals=unique(ID);
    ind_col=interp1(col,size(col,1)/length(individuals):size(col,1)/length(individuals):size(col,1),[],'extrap');
    spreadIDrange=[-0.1:0.01:0.1];
    spreadID=interp1(spreadIDrange,length(spreadIDrange)/length(individuals):length(spreadIDrange)/length(individuals):length(spreadIDrange),[],'extrap');

    subplot(2,4,1);hold on;
    b=boxplot(positions(:,1),'color',boxCol);set(b,{'linew'},{2});hold on;
    for u=1:length(individuals)
        temp=find(ID==individuals(u));
        n_u=length(temp);
        plot(spreadID(u)+ones(n_u,1),positions(temp,1),'.','color',ind_col(u,:),'MarkerSize',16);hold on;
    end
    % plot(0.05*randn(n,1)+ones(n,1),perc_half(:,2),'.r','MarkerSize',16);
    plot([0.5 1.5],[0 0 ],'k--');
    title('median pos');ylim([-150 150]);xlim([xmin xmax]);
    ylabel('mm');

    subplot(2,4,2);hold on;
    b=boxplot(speeds(:,1),'color',boxCol);set(b,{'linew'},{2});hold on;
    for u=1:length(individuals)
        temp=find(ID==individuals(u));
        n_u=length(temp);
        plot(spreadID(u)+ones(n_u,1),speeds(temp,1),'.','color',ind_col(u,:),'MarkerSize',16);hold on;
    end
    title('avg speed');ylim([0 2500]);xlim([xmin xmax]);
    ylabel('mm/s');

    tunnelCenter=[[tunnel(1,1), tunnel(1,2)-0.5*(tunnel(1,2) - tunnel(3,2))];[tunnel(2,1), tunnel(1,2)-0.5*(tunnel(1,2) - tunnel(3,2))]];

    subplot(2,4,3);hold on;
    b=boxplot(positions(:,2),'color',boxCol);set(b,{'linew'},{2});hold on;
    for u=1:length(individuals)
        temp=find(ID==individuals(u));
        n_u=length(temp);
        plot(spreadID(u)+ones(n_u,1),positions(temp,2),'.','color',ind_col(u,:),'MarkerSize',16);hold on;
    end
    title('pos variation');ylim([0 100]);xlim([xmin xmax]);
    ylabel('std (mm)');

    subplot(2,4,4);hold on;
    b=boxplot(tortuosity,'color',boxCol);set(b,{'linew'},{2});hold on;
    for u=1:length(individuals)
        temp=find(ID==individuals(u));
        n_u=length(temp);
        plot(spreadID(u)+ones(n_u,1),tortuosity(temp),'.','color',ind_col(u,:),'MarkerSize',16);hold on;
    end
    title('avg tortuosity');ylim([1 1.5]);xlim([xmin xmax]);
    ylabel('length/distance');

    subplot(2,4,5);hold on;
    b=boxplot(positions(:,3),'color',boxCol);set(b,{'linew'},{2});hold on;
    for u=1:length(individuals)
        temp=find(ID==individuals(u));
        n_u=length(temp);
        plot(spreadID(u)+ones(n_u,1),positions(temp,3),'.','color',ind_col(u,:),'MarkerSize',16);hold on;
    end
    title('crossover index');ylim([-150 150]);xlim([xmin xmax]);
    ylabel('std (mm/s)');

    subplot(2,4,6);hold on;
    b=boxplot(positions(:,5),'color',boxCol);set(b,{'linew'},{2});hold on;
    for u=1:length(individuals)
        temp=find(ID==individuals(u));
        n_u=length(temp);
        plot(spreadID(u)+ones(n_u,1),positions(temp,5),'.','color',ind_col(u,:),'MarkerSize',16);hold on;
    end
    % plot(0.05*randn(n,1)+ones(n,1),perc_half(:,2),'.r','MarkerSize',16);
    b=boxplot(positions(:,5),'color',boxCol);set(b,{'linew'},{2});hold on;
    title('prop lateral mov');ylim([0 2]);xlim([xmin xmax]);
    ylabel('y/x travel');

    subplot(2,4,7);hold on;
    b=boxplot(positions(:,4),'color',boxCol);set(b,{'linew'},{2});hold on;
    for u=1:length(individuals)
        temp=find(ID==individuals(u));
        n_u=length(temp);
        plot(spreadID(u)+ones(n_u,1),positions(temp,4),'.','color',ind_col(u,:),'MarkerSize',16);hold on;
    end
    plot([0.5 1.5],[0 0 ],'k--');
    title('interp med pos');ylim([-150 150]);xlim([xmin xmax]);
    ylabel('mm');


    subplot(2,4,8);hold on
    b=boxplot(speeds(:,2),'color',boxCol);set(b,{'linew'},{2});hold on;
    for u=1:length(individuals)
        temp=find(ID==individuals(u));
        n_u=length(temp);
        plot(spreadID(u)+ones(n_u,1),speeds(temp,2),'.','color',ind_col(u,:),'MarkerSize',16);hold on;
    end
    title('max 15pc speed');ylim([0 3000]);xlim([xmin xmax]);
    ylabel('mm/s');

end
end