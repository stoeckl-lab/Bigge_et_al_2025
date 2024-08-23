clear
close all
%generate data table for statistical analysis of all optic flow and
%contrast data

%% loading files
% startDir='Y:\Anna&Ronja\Dorsal_Ventral_Imaging\ContrastAnalysis\All';%old
% startDir='F:\Anna Backup\AG Stoeckl\Anna&Ronja\Dorsal_Ventral_Imaging\ResultsNew\Natural';%new
% startDir='F:\Anna Backup\AG Stoeckl\Anna&Ronja\Dorsal_Ventral_Imaging\ResultsNew\Natural\of';%new

startDir='F:\Anna Backup\AG Stoeckl\Anna&Ronja\Dorsal_Ventral_Imaging\ResultsNew\Tunnel\contrast';%new
startDir='F:\Anna Backup\AG Stoeckl\Anna&Ronja\Dorsal_Ventral_Imaging\ResultsNew\Tunnel\of';%new

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

    %check habitat conditions
    if (contains(filenames{i},'Open') && ~ contains(filenames{i},'Semi')) || contains(filenames{i},'Flower')
        condition{i}='open';
    elseif contains(filenames{i},'Forest')
        condition{i}='closed';
    elseif contains(filenames{i},'Semi')
        condition{i}='semi';
    end

    %save scene type
    ind1=find(filenames{i}=='_',1,'first');
    scene{i}=filenames{i}(1:ind1-1);

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
t=table(condition',scene',magnitude(:,1),magnitude(:,2),magnitude(:,3),magnitude(:,4),'VariableNames',{'condition','scene','dorsal','ventral','left','right'});

if contains(filenames{i},'contrast')
    writetable(t,'contrast_data_all.xls');
else
    writetable(t,'of_data_all.xls');
end

%make long format with quadrants as factors
t1=table(t{:,1},t{:,2},t{:,3},repmat({'dorsal'},size(t{:,1},1),1),'VariableNames',{'condition','scene','magnitude','quadrant'});
t2=table(t{:,1},t{:,2},t{:,4},repmat({'ventral'},size(t{:,1},1),1),'VariableNames',{'condition','scene','magnitude','quadrant'});
t3=table(t{:,1},t{:,2},t{:,5},repmat({'left'},size(t{:,1},1),1),'VariableNames',{'condition','scene','magnitude','quadrant'});
t4=table(t{:,1},t{:,2},t{:,6},repmat({'right'},size(t{:,1},1),1),'VariableNames',{'condition','scene','magnitude','quadrant'});

tLong=[t1;t2;t3;t4];

if contains(filenames{i},'contrast')
    writetable(tLong,'contrast_data_all_long.xls');
else
    writetable(tLong,'of_data_all_long.xls');
end

%save long format with lateral (average left and right)
t3=table(t{:,1},t{:,2},0.5*(t{:,5}+t{:,6}),repmat({'lateral'},size(t{:,1},1),1),'VariableNames',{'condition','scene','magnitude','quadrant'});

tLongLateral=[t1;t2;t3];
if contains(filenames{i},'contrast')
    writetable(tLongLateral,'contrast_data_all_long_lateral.xls');
else
    writetable(tLongLateral,'of_data_all_long_lateral.xls');
end
