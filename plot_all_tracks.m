function plot_all_tracks 
%load all digitised flight tracks ("auto_data...mat") in custom selection
% and plot them superimposed, with median and interquartile ranges of the
% population


%load filename with xypts extension
files=uigetfile('multiselect','on','*.mat');
if iscell(files)==0
    temp=files;
    clear files
    files=cell(1,1); files{1}=temp;
end

%if I load comp_measures file or if I load all_tracks file
if isempty(strfind(files{1},'comp_measures'))==0 || isempty(strfind(files{1},'all_tracks'))==0
    load(files{1})
    if exist('allTracks')
        clear files
        files=allTracks;
    end
end


%plot color
plotColor=[0.5 0.5 0.5];

allTracksInterp=nan(length(files),1000);
allSpeedsInterp=nan(length(files),1000);
%plot simple lines
f1=figure;hold on;

for i=1:length(files)
    if exist('allTracks')==0
        load(files{i});
    end
    
    %     %make all start from left, i.e. flip moths that cross tunnel from right
    %         if strfind(files{i},'data_R')
    %             moth(:,1)=1000- moth(:,1);
    %         end
    
    %generate matrix to plot heatmap
    %some traces have the upper tunnel border up, and the lower one down.
    %have to swap them!
    if exist('tunnel')
        if tunnel(3,2)<tunnel(1,2)
            temp=tunnel(1:4,2);tunnel(3:4,2)=temp(1:2);tunnel(1:2,2)=temp(3:4);
        end
    else tunnel=[0 0;0 0;0 0;0 0];%in this case, load compe_measures, where traces are already aligned to tunnel proportions
    end
    
    if exist('moth') %this is for loading tunnel_track data
        if length(moth(:,2))==length(unique(moth(:,2))) && length(moth(:,1))==length(unique(moth(:,1)))
            %simple line plot
            moth(:,2)=moth(:,2)-nanmean(tunnel(1:4,2));
            plot(moth(:,1),moth(:,2),'-','color',plotColor);hold on;%legend({'flower', 'moth'})
            set(gca,'PlotBoxAspectRatio',[3.33 1 1])
            
            %collect all tracks interpolated to make average and range
            allTracksInterp(i,:)=interp1(moth(~isnan(moth(:,1)),1),moth(~isnan(moth(:,1)),2),[1:1:1000]);
        end
    elseif exist('allTracks') %this is for loading comp_measures
        
        if length(allTracks{i}(:,2))==length(unique(allTracks{i}(:,2))) && length(allTracks{i}(:,1))==length(unique(allTracks{i}(:,1)))
            plot(allTracks{i}(:,1),allTracks{i}(:,2),'-','color',plotColor);hold on;%legend({'flower', 'moth'})
            set(gca,'PlotBoxAspectRatio',[3.33 1 1])
            
            %collect all tracks interpolated to make average and range
            
            allTracksInterp(i,:)=interp1(allTracks{i}(:,1),allTracks{i}(:,2),[1:1:1000]);
            allSpeedsInterp(i,:)=interp1(all_pos{i}(:,1),all_speed{i}(:,1),[1:1:1000]);
        end
    end
    
end

plot([0,1000], [0 0],'k--');
xlim([0 1000])
ylim([-150 150])
title(['n=',num2str(length(files))]);
box on

% hold on
medianData=nanmedian(allTracksInterp,1);
rangeDataUp=quantile(allTracksInterp,.75)-quantile(allTracksInterp,.5);
rangeDataDown=quantile(allTracksInterp,.5)-quantile(allTracksInterp,.25);

% rangeData=nanstd(allTracksInterp,1,1)/(sqrt(size(allTracksInterp,1)));

indDel=find(abs(medianData)>2*quantile(abs(medianData),.75));
medianData(indDel)=nan;rangeData(indDel)=nan;

shadedErrorBar_anna([1:1:1000],smooth(medianData,5),[smooth(rangeDataUp,5)';smooth(rangeDataDown,5)'],'r')
tunnelCenter=[[tunnel(1,1), tunnel(1,2)-0.5*(tunnel(1,2) - tunnel(3,2))];[tunnel(2,1), tunnel(1,2)-0.5*(tunnel(1,2) - tunnel(3,2))]];

print -f1 -dpdf -r300 -painters all_tracks.eps

end