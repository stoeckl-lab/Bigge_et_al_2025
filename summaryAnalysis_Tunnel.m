clear
close all
%generate data table for statistical analysis of all optic flow and
%contrast data

%% loading files

startDir='F:\Anna Backup\AG Stoeckl\Anna&Ronja\Dorsal_Ventral_Imaging\ResultsNew\Tunnel\contrast';%new
% startDir='F:\Anna Backup\AG Stoeckl\Anna&Ronja\Dorsal_Ventral_Imaging\ResultsNew\Tunnel\of';%new

cd(startDir)

filenames=uigetfile('*mat','MultiSelect','on');
if ischar(filenames)
    temp=filenames;
    filenames=cell(1,1);
    filenames{1}=temp;
end

condition=cell(size(filenames));
magnitude=nan(length(filenames),4);
scene=cell(size(filenames));

for i=1:length(filenames)
    load(filenames{i})

    %check tunnel conditions
    ind1=find(filenames{i}=='_');
    ind2=find(filenames{i}=='.',1,'last');
    condition{i}=filenames{i}(ind1+1:ind2-1);
    %remove "contrast" from condition
    indDel=strfind(condition{i},'contrast');
    if ~isnan(indDel)
    condition{i}(indDel:end)=[];
    end

    if exist("cdata")
        if exist("cdata.MeanMagnitude")
            magnitude(i,:)=cdata.MeanMagnitude;
        else
            %calculate conrtast in the quartiles
            lowerTriag=tril(cdata.MedianMag_all);
            upperTriag=triu(cdata.MedianMag_all);

            leftQ=tril(flipud(lowerTriag));
            rightQ=(triu(flipud(upperTriag)));

            ventralQ=flipud(triu(flipud(lowerTriag)));
            dorsalQ=flipud(tril(flipud(upperTriag)));
            magnitude(i,:)=[nanmean(dorsalQ(:)) nanmean(ventralQ(:)) nanmean(leftQ(:)) nanmean(rightQ(:))];
        end
    else
        magnitude(i,:)=data.MeanMagnitude;
    end
end


%% save data in table
%collect data in table
t=table(condition',magnitude(:,1),magnitude(:,2),magnitude(:,3),magnitude(:,4),'VariableNames',{'condition','dorsal','ventral','left','right'});

if contains(filenames{i},'contrast')
    writetable(t,'contrast_data_all.xls');
else
    writetable(t,'of_data_all.xls');
end

%make long format with quadrants as factors
t1=table(t{:,1},t{:,2},repmat({'dorsal'},size(t{:,1},1),1),'VariableNames',{'condition','magnitude','quadrant'});
t2=table(t{:,1},t{:,3},repmat({'ventral'},size(t{:,1},1),1),'VariableNames',{'condition','magnitude','quadrant'});
t3=table(t{:,1},t{:,4},repmat({'left'},size(t{:,1},1),1),'VariableNames',{'condition','magnitude','quadrant'});
t4=table(t{:,1},t{:,5},repmat({'right'},size(t{:,1},1),1),'VariableNames',{'condition','magnitude','quadrant'});

tLong=[t1;t2;t3;t4];

if contains(filenames{i},'contrast')
    writetable(tLong,'contrast_data_all_long.xls');
else
    writetable(tLong,'of_data_all_long.xls');
end

%save long format with lateral (average left and right)
t3=table(t{:,1},0.5*(t{:,4}+t{:,5}),repmat({'lateral'},size(t{:,1},1),1),'VariableNames',{'condition','magnitude','quadrant'});

tLongLateral=[t1;t2;t3];
if contains(filenames{i},'contrast')
    writetable(tLongLateral,'contrast_data_all_long_lateral.xls');
else
    writetable(tLongLateral,'of_data_all_long_lateral.xls');
end
