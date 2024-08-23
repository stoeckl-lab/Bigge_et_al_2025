clear
close all
%plot selected optic flow or contrast data, average when multiple scenes
%selected

%% loading files
% startDir='Y:\Anna&Ronja\Dorsal_Ventral_Imaging\ContrastAnalysis\All';%old
% startDir='F:\Anna Backup\AG Stoeckl\Anna&Ronja\Dorsal_Ventral_Imaging\ResultsNew\Natural\contrast';%new
% startDir='F:\Anna Backup\AG Stoeckl\Anna&Ronja\Dorsal_Ventral_Imaging\ResultsNew\Natural\of';%new

% startDir='F:\Anna Backup\AG Stoeckl\Anna&Ronja\Dorsal_Ventral_Imaging\ResultsNew\Tunnel\contrast';%new
startDir='F:\Anna Backup\AG Stoeckl\Anna&Ronja\Dorsal_Ventral_Imaging\ResultsNew\Tunnel\of';%new

cd(startDir)

filenames=uigetfile('*mat','MultiSelect','on');
if ischar(filenames)
    temp=filenames;
    filenames=cell(1,1);
    filenames{1}=temp;
end

allData=nan(1200,1200,length(filenames));
contrastDataNorm=nan(1200,1200,length(filenames));

for i=1:length(filenames)
    load(filenames{i})
    if exist("cdata")
        allData(:,:,i)=cdata.MedianMag_all;
    else
        allData(:,:,i)=data.MedianMag_all;
    end

    %for the supplementary data in the paper
    tempData=allData(:,:,i);
    contrastDataNorm(:,:,i)=allData(:,:,i)./quantile(tempData(:),0.995);
end

%% Plotting all in one figure
if contains(filenames{1},'contrast')
    isContrast=true;
    %make a "cool" colormap
    % cols=[0 0 0;68 85 196;69 137 252;42 239 160;165 252 69;255 255 255]/255;
    cols=[0 0 0;68 85 196;69 137 252; 39 170 225;39 255 255;180 255 255;255 255 255]/255;
    colmap=interp1(ceil((0:size(cols,1)-1)*(256/(size(cols,1)-1))),cols,1:256);
    zscale=[0 0.003];
else
    isContrast=false;
    colmap=colormap('hot');
    zscale=[0 3];
end

f1= figure;
meanData=nanmean(allData,3);
% meanData=nanmean(contrastDataNorm,3);
%normalise to one (though take care of very few but high outliers by using
%quantile 0.995
% meanData=meanData/quantile(meanData(:),0.99);
% imagesc(meanData,[0 1])
imagesc(meanData,zscale)

colormap(colmap)
colorbar

hold on
% plot([0 1200],[600 600], 'w', 'LineWidth', 1.5);
axis equal
axis off

% save figure
%hack figure names for compound figures
% filenames{1}='lateral.mat';

if isContrast
    print('-f1','-dpdf', '-r300', '-painters', '-bestfit', [filenames{1}(1:end-4),'_contrast_heatmap.eps'])
else
    print('-f1','-dpdf', '-r300', '-bestfit', [filenames{1}(1:end-4),'_heatmap.eps'])
end

