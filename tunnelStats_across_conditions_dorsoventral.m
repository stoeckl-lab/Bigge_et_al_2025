%plot figures and statistical results across tunnel conditions

clear all;close all


%% select input directories
cd 'flightData'

%select a variable number of .mat files in the right order
% if exist('inputDirs')==0
% count=0;
% uiind=1;
% inputDirs=cell(length(dir)-2,1);%give it maximum length
% while uiind==1
% tempDirs=uigetdir;
% count=count+1;
% inputDirs{count}=tempDirs;
% if tempDirs==0
%     uiind=0;
% end
% end
% inputDirs(count:end)=[];%delete empty entries
% end

% resDirName='Fig2D'; inputDirs={'none', 'rightOF', 'leftOF','lateralOF','ventralOF','dorsalOF'};
% resDirName='Fig2G'; inputDirs={'ventralSwitch_rl', 'ventralSwitch_lr','dorsalSwitch_rl','dorsalSwitch_lr','dorsal_switchgratingLR','dorsal_switchgratingRL'};
% resDirName='Fig2H'; inputDirs={'ventral_halfOF_l', 'ventral_halfOF_r','dorsal_halfLong_l','dorsal_halfOF_r','dorsal_halfOF_r_long_l','dorsal_halfLong_r','dorsal_halfOF_l','dorsal_halfOF_l_long_r'};
% resDirName='Fig2I'; inputDirs={'dorsalOF_3cm', 'dorsalOF', 'dorsalOF_12cm'};
% resDirName='Fig3'; inputDirs={'leftOF', 'dorsalSwitch_rl','dorsalSwitch_rl_leftOF','rightOF','dorsalSwitch_lr','dorsalSwitch_lr_rightOF'};
% resDirName='Fig4C'; inputDirs={'rightOF', 'dorsal_halfOF_l','dorsal_halfOF_l_rightOF','leftOF','dorsal_halfOF_r','dorsal_halfOF_r_leftOF'};
% resDirName='Fig2D'; inputDirs={'ventralOF2','dorsalOF2','dorsalOF_ventralOF'};
% resDirName='FigS2'; inputDirs={'lateralChecker', 'dorsalSwitch_rl','lateralChecker_switch_rl','dorsalSwitch_lr','lateralChecker_switch_lr'};
resDirName='FigS3C'; inputDirs={'none2', 'ventral_long','ventralOF2','dorsal_long','dorsalOF2','dorsalOF_ventralOF'};

%% initialise variables
n=length(inputDirs);

all_avg_speed=[];
all_max_speed=[];
% all_var_speed=[];
all_tortuosity=[];
all_position=[];
all_var_position=[];
all_cross_position=[];
all_lateral_position=[];
all_mean_areas=[];
all_var_areas=[];
% all_position_1stH=[];
% all_cross_position_1stH=[];
all_tracks=[];

groups=[];

%% loop through directories and extract data

for i=1:length(inputDirs)
%     cd(inputDirs{i})
%     load comp_measures.mat 
%load .mat files with data
    load([inputDirs{i},'.mat']);

    %extract speed
    speedPercentage=0.1; %what percentage of all speeds as maximum speed
    all_avg_speed=[all_avg_speed; speeds(:,1)];
    for t=1:length(all_speed)
        if isempty(nanmax(all_speed{t}))
            temp_max_speed(t)=nan; 
        else
        temp_max_speed(t)=quantile(all_speed{t},1-speedPercentage);%compute the max. x% speed
        end
    end
    all_max_speed=[all_max_speed; temp_max_speed'];%speeds(:,2)
    clear temp_max_speed;
%     all_var_speed=[all_var_speed; speeds(:,3)];
    
%     rel_histSpeedData=histSpeedS.counts./repmat(nansum(histSpeedS.counts,2),1,size(histSpeedS.counts,2));
%     all_histSpeed_mean(i,:)=nanmean(rel_histSpeedData,1);
%     all_histSpeed_conf(i,:)=nanstd(rel_histSpeedData,1,1)/(sqrt(size(histSpeedS.counts,1))-1);
%     all_all_histSpeedEdges=histSpeedS.edges(1:end-1)+0.5*nanmean(diff(histSpeedS.edges));

    %tortuosity
    all_tortuosity=[all_tortuosity; tortuosity];
    
    %position
    all_position=[all_position; positions(:,1)];
    all_var_position=[all_var_position; positions(:,2)];
    all_cross_position=[all_cross_position; positions(:,3)];
    all_lateral_position=[all_lateral_position; positions(:,5)];
%     all_position_1stH=[all_position_1stH; positions(:,6)];
%     all_cross_position_1stH=[all_cross_position_1stH; positions(:,7)];
if exist('allTracks') %some of the older analysis files don't have tracks
    all_tracks=[all_tracks; allTracks];%collect all individual flight tracks
end
    %area, if it exists
    if exist('all_areas')
        if iscell(all_areas)
        tempMean=nan(length(all_areas),3);tempVar=nan(length(all_areas),3);
        for u=1:length(all_areas)
            if ~isempty(all_areas{u})
            tempMean(u,:)=nanmean(all_areas{u},1);
            tempVar(u,:)=nanstd(all_areas{u},1,1);
            else
            tempMean(u,:)=[nan nan nan];
            tempVar(u,:)=[nan nan nan];
            end
        end        
        all_mean_areas=[all_mean_areas; tempMean];
        all_var_areas=[all_var_areas; tempVar];
        end
    end
    
    groups=[groups;i*ones(size(speeds(:,1)))];
    
    if isempty(strfind(inputDirs{i},'\'))
        label{i}=inputDirs{i};
    else
    indStr=strfind(inputDirs{i},'\');
    label{i}=inputDirs{i}(indStr(end)+1:end);
    end
%     cd ..
end


%% plot data

%cd into results dir to save plots and stats there
cd D:\Nextcloud\Home\Behaviour\Tunnel\2018-2020_dorsal_ventral_Maximilian_Rebecca\Results_DV_2024
if ~isfolder(resDirName)
    mkdir(resDirName)
end
cd(resDirName)

%set plot parameters
xmin=0.8;xmax=1.2;spread=0.05;
%make dot colour rainbow depending on number of entries
dotCol=[0 0 0];
boxCol=[1 1 1];
violinCol=[168 211 65]/255;
%rotate labels to be angles at 0 to not scale the figures

%adjust the width of the figures by the number of groups
width=90*length(unique(groups));

%POSITION
plot_label='median pos';
%this is what we used for dorsal ventral CB paper
bw=range(all_position(:))/25;%bandwidth for violin plots
%HOW TO USE DEFAULT?

%initialise stats file
diary stats.txt

f1=figure('Position',[500 500 width 400]);hold on;

for u=1:length(unique(groups))
    n=sum(groups==u);
v=violin(all_position(groups==u),'x',u,'facecolor',violinCol,'edgecolor',[],'medc','w','facealpha',0.8,'bw',bw );
plot(spread*randn(n,1)+u*ones(n,1),all_position(groups==u),'.','MarkerSize',12,'color',dotCol);hold on
end
b=boxplot(all_position,groups,'color',boxCol,'PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol',' ','LabelOrientation','horizontal');
% set(b,{'linew'},{2})
plot([0.5 length(label)+0.5],[0 0],'k--');
title(plot_label);ylim([-150 150]);%xlim([0.5 1.5]);
% set(gca,'xticklabel',label);
ylabel([plot_label '(mm)'])


group_stats(all_position, groups,plot_label)

%pairwise comparisons using vartest
disp('median position vartest');
ind_all=nchoosek(unique(groups),2);
for i=1:size(ind_all,1)
    ind=[find(groups==ind_all(i,1)); find(groups==ind_all(i,2))];
[p,stats] = vartestn(all_position(ind),groups(ind),'TestType','BrownForsythe','Display','off');
table(i,1:2)=ind_all(i,:);
table(i,3:4)=[p stats.fstat];
end
disp(table);

%test each condition against 0 mean
disp('median position 0 mean');
for i=1:max(groups)
    ind=find(groups==i);
[p,stats] = signrank(all_position(ind));
table2(i,1)=i;
table2(i,2)=[stats];
table2(i,3)=[p];
end
disp(table2);


%AVG SPEED
plot_label='avg speed';
bw=range(all_avg_speed(:))/25;%bandwidth for violin plots

f2=figure('Position',[500 500 width 400]);hold on;
for u=1:length(unique(groups))
    n=sum(groups==u);
v=violin(all_avg_speed(groups==u),'x',u,'facecolor',violinCol,'edgecolor',[],'medc','w','facealpha',0.8,'bw',bw);
plot(spread*randn(n,1)+u*ones(n,1),all_avg_speed(groups==u),'.','MarkerSize',12,'color',dotCol);hold on
end
b=boxplot(all_avg_speed,groups,'color',boxCol,'PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol',' ','LabelOrientation','horizontal');
% b=boxplot(all_avg_speed,groups,'color',boxCol);set(b,{'linew'},{2})
title('avg speed');ylim([200 2200]);%xlim([0.5 1.5]);
ylabel([plot_label ' (mm/s)']);
% set(gca,'xticklabel',label);

group_stats(all_avg_speed, groups,plot_label)


%POSITION VARIATION
plot_label='var. in pos.';

f3=figure('Position',[500 500 width 400]);hold on;
for u=1:length(unique(groups))
    n=sum(groups==u);
v=violin(all_var_position(groups==u),'x',u,'facecolor',violinCol,'edgecolor',[],'medc','w','facealpha',0.8);
plot(spread*randn(n,1)+u*ones(n,1),all_var_position(groups==u),'.','MarkerSize',12,'color',dotCol);hold on
end
b=boxplot(all_var_position,groups,'color',boxCol,'PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol',' ','LabelOrientation','horizontal');
% b=boxplot(all_var_position,groups,'color',boxCol);set(b,{'linew'},{2})
title(plot_label);ylim([0 100]);%xlim([0.5 1.5]);
set(gca,'xticklabel',label);
ylabel([plot_label ' (mm)'])

group_stats(all_var_position, groups,plot_label)

% 
%CROSS POSITION
plot_label='cross pos index';
bw=range(all_cross_position(:))/25;%bandwidth for violin plots

f4=figure('Position',[500 500 width 400]);hold on;
for u=1:length(unique(groups))
    n=sum(groups==u);
    v=violin(all_cross_position(groups==u),'x',u,'facecolor',violinCol,'edgecolor',[],'medc','w','facealpha',0.8,'bw',bw);
plot(spread*randn(n,1)+u*ones(n,1),all_cross_position(groups==u),'.','MarkerSize',12,'color',dotCol);hold on
end
b=boxplot(all_cross_position,groups,'color',boxCol,'PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol',' ','LabelOrientation','horizontal');
% b=boxplot(all_cross_position,groups,'color',boxCol);set(b,{'linew'},{2})
plot([0.5 length(label)+0.5],[0 0],'k--');
title(plot_label);ylim([-200 200]);%xlim([0.5 1.5]);
% set(gca,'xticklabel',label);
ylabel([plot_label ' (mm)'])

group_stats(all_cross_position, groups,plot_label)


%LATERAL MOVEMENT
plot_label='prop lateral movement';
bw=range(all_lateral_position(:))/25;%bandwidth for violin plots

f5=figure('Position',[500 500 width 400]);hold on;
for u=1:length(unique(groups))
    n=sum(groups==u);
v=violin(all_lateral_position(groups==u),'x',u,'facecolor',violinCol,'edgecolor',[],'medc','w','facealpha',0.8,'bw',bw);
plot(spread*randn(n,1)+u*ones(n,1),all_lateral_position(groups==u),'.','MarkerSize',12,'color',dotCol);hold on
end
b=boxplot(all_lateral_position,groups,'color',boxCol,'PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol',' ','LabelOrientation','horizontal');
% b=boxplot(all_lateral_position,groups,'color',boxCol);set(b,{'linew'},{2})
title(plot_label);ylim([0 1.5]);%xlim([0.5 1.5]);
% set(gca,'xticklabel',label);
ylabel('y/x movement')

group_stats(all_lateral_position, groups,plot_label)


%MAX SPEED
plot_label=sprintf('max %u %% speed',speedPercentage*100);
f6=figure('Position',[500 500 width 400]);hold on;
for u=1:length(unique(groups))
    n=sum(groups==u);
v=violin(all_max_speed(groups==u),'x',u,'facecolor',violinCol,'edgecolor',[],'medc','w','facealpha',0.8);
plot(spread*randn(n,1)+u*ones(n,1),all_max_speed(groups==u),'.','MarkerSize',12,'color',dotCol);hold on
end
b=boxplot(all_max_speed,groups,'color',boxCol,'PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol',' ','LabelOrientation','horizontal');
title(plot_label);
ylim([0 2500]);%xlim([0.5 1.5]);
ylabel('mm/s');
% set(gca,'xticklabel',label);

group_stats(all_max_speed, groups,plot_label)

% %% save data for the spatial response tunnel paper
% data=[all_position all_var_position all_cross_position all_lateral_position all_tortuosity all_avg_speed all_max_speed];
% dataLabels={'median position','in flight variation','cross position','rel lateral movement','tortuosity','avg speed','max speed'};
% conditions=label;
% save('all_raw_data.mat','data','dataLabels','groups','conditions','all_tracks');

%% for the spatial response paper

%save data from both the symmetric and asymmetric condition to load here
%again
% data=[all_position all_var_position all_cross_position all_lateral_position all_tortuosity all_avg_speed all_max_speed];
% dataLabels={'median position','in flight variation','cross position','rel lateral movement','tortuosity','avg speed','max speed'};
% conditions=label;
% save('allData.mat','data','dataLabels','groups','conditions','all_avg_speed','all_max_speed','all_tortuosity','all_position','all_var_position','all_cross_position','all_lateral_position');

% %plot flight speed vs lateral movement
% figure;hold on;
% plot(log(all_avg_speed),log(all_lateral_position),'.','markersize',12);hold on;
% axis square;xlabel('log speed');ylabel('log lateral movement')
% %remove nan values
% x1=log(all_avg_speed);
% x2=log(all_lateral_position);
% indnan=isnan(x1)+isnan(x2);
% [r1,m1,b1]=regression(x1(indnan==0)',x2(indnan==0)');
% [rho1,p1]=corrcoef(x1(indnan==0),x2(indnan==0));
% plot(x1,polyval([m1,b1],x1),'k','LineWidth',2);
% title(sprintf('r= %1.2f, p=%1.3f, m=%1.3f, b=%1.3f',r1,p1(2),m1,b1));


%%
% % POSITION in 1st half
% plot_label='median pos 1st half';
% f7=figure('Position',[500 500 width 400]);hold on;
% for u=1:length(unique(groups))
%     n=sum(groups==u);
% v=violin(all_position_1stH(groups==u),'x',u,'facecolor',violinCol,'medc','w','facealpha',0.8);
% plot(spread*randn(n,1)+u*ones(n,1),all_position_1stH(groups==u),'.','MarkerSize',12,'color',dotCol);hold on
% end
% b=boxplot(all_position_1stH,groups,'color',boxCol,'PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol',' ','LabelOrientation','horizontal');
% % set(b,{'linew'},{2})
% plot([0.5 length(label)+0.5],[0 0],'k--');
% title(plot_label);ylim([-150 150]);%xlim([0.5 1.5]);
% set(gca,'xticklabel',label);
% ylabel([plot_label ' (mm)'])
% 
% group_stats(all_position_1stH, groups,plot_label)

% 
% %test each condition against 0 mean
% disp('median position 1st half 0 mean');
% for i=1:max(groups)
%     ind=find(groups==i);
% [p,stats] = signrank(all_position_1stH(ind));
% table2(i,1)=i;
% table2(i,2)=[stats];
% table2(i,3)=[p];
% end
% disp(table2);


%CROSS POSITION 1st half
% plot_label='cross position index 1stH';
% f8=figure('Position',[500 500 width 400]);hold on;
% for u=1:length(unique(groups))
%     n=sum(groups==u);
%     v=violin(all_cross_position_1stH(groups==u),'x',u,'facecolor',violinCol,'medc','w','facealpha',0.8);
% plot(spread*randn(n,1)+u*ones(n,1),all_cross_position_1stH(groups==u),'.','MarkerSize',12,'color',dotCol);hold on
% end
% b=boxplot(all_cross_position_1stH,groups,'color',boxCol,'PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol',' ','LabelOrientation','horizontal');
% % b=boxplot(all_cross_position,groups,'color',boxCol);set(b,{'linew'},{2})
% plot([0.5 length(label)+0.5],[0 0],'k--');
% title(plot_label);ylim([-200 200]);%xlim([0.5 1.5]);
% set(gca,'xticklabel',label);
% ylabel([plot_label ' (mm)'])
% 
% group_stats(all_cross_position_1stH, groups,plot_label)

%% PLOT AREA OF TRACKING
% 
%check if all data have the area parameter, if not, skip plot
try grpstats(all_mean_areas,groups,'sum')
% MEAN AREA
plot_label=sprintf('area',speedPercentage*100);

f9=figure('Position',[500 500 width 400]);hold on;
for u=1:length(unique(groups))
    n=sum(groups==u);
v=violin(all_mean_areas(groups==u,1),'x',u,'facecolor',violinCol,'medc','w','facealpha',0.8);
plot(spread*randn(n,1)+u*ones(n,1),all_mean_areas(groups==u,1),'.','MarkerSize',12,'color',dotCol);hold on
end
b=boxplot(all_mean_areas(:,1),groups,'color',boxCol,'PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol',' ','LabelOrientation','horizontal');
title(plot_label);
% ylim([0 2500]);%xlim([0.5 1.5]);
ylabel([plot_label '(mm^2)']);
% set(gca,'xticklabel',label);

group_stats(all_mean_areas, groups,plot_label)


% VARIATION IN AREA
plot_label=sprintf('variation in area',speedPercentage*100);

f10=figure('Position',[500 500 width 400]);hold on;
for u=1:length(unique(groups))
    n=sum(groups==u);
v=violin(all_var_areas(groups==u,1),'x',u,'facecolor',violinCol,'medc','w','facealpha',0.8);
plot(spread*randn(n,1)+u*ones(n,1),all_var_areas(groups==u,1),'.','MarkerSize',12,'color',dotCol);hold on
end
b=boxplot(all_var_areas(:,1),groups,'color',boxCol,'PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol',' ','LabelOrientation','horizontal');
title(plot_label);
% ylim([0 2500]);%xlim([0.5 1.5]);
ylabel(plot_label);
set(gca,'xticklabel',label);

group_stats(all_var_areas, groups,plot_label)


% %area histogram 
% area_sizes=[50:5:300];
% f11=figure('Position',[500 500 width 400]);hold on;
% for u=1:length(unique(groups))
% histogram(all_var_areas(groups==u,1),area_sizes);hold on
% end
catch
end

% end



%stop writing stats into file
diary off
%% save all data

%% save plots
% cd Results
% %can be tricky when file size too large (though not size, not sure what it
% %is?)
% %read more here: https://de.mathworks.com/matlabcentral/answers/92521-why-does-matlab-not-export-eps-files-properly
set(f1,'PaperPositionMode','auto');print -f1 -dpdf -r300 -painters medPosition.eps
set(f2,'PaperPositionMode','auto');print -f2 -dpdf -r300 -painters avgSpeed.eps
set(f3,'PaperPositionMode','auto');print -f3 -dpdf -r300 -painters varPosition.eps
set(f4,'PaperPositionMode','auto');print -f4 -dpdf -r300 -painters crossPosition.eps
set(f5,'PaperPositionMode','auto');print -f5 -dpdf -r300 -painters lateralMov.eps
set(f6,'PaperPositionMode','auto');print -f6 -dpdf -r300 -painters maxSpeed.eps
% set(f7,'PaperPositionMode','auto');print -f7 -dpdf -r300 -painters pos1stH.eps
% set(f8,'PaperPositionMode','auto');print -f8 -dpdf -r300 -painters crossPos1stH.eps

if exist('f9')
set(f9,'PaperPositionMode','auto');print -f7 -dpdf -r300 -painters area.eps
set(f10,'PaperPositionMode','auto');print -f8 -dpdf -r300 -painters varArea.eps
end

% dir=uigetdir;
% cd(dir);
% %% 

%% inbuilt functions

function group_stats(data, groups,label)

%test normality
residuals=[];
[p,tbl,stats] = anova1(data(:,1),groups,'off');
c = multcompare(stats,'display','off');
for u=1:length(unique(groups))
residuals=[residuals;data(groups==u,1)-stats.means(u)];
end
%choose right test
if lillietest(residuals)==0
disp('ANOVA');
disp(label);
disp(tbl)
disp(c(:,[1,2,6]));
else
disp('KRUSKAL WALLIS');
[p,tbl,stats] = kruskalwallis(data(:,1),groups,'off');
c = multcompare(stats,'display','off');
disp(label);
disp(tbl)
disp(c(:,[1,2,6]));
end% 
end
