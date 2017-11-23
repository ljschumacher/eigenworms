% plots the distribution of eigenworm projections for in and out-of-cluster
% worms
close all
clear

% issues/to-do:
% - filter red channel data further? currently we are assuming that those
% with feature calculations are filtered fairly well

% figure export options
exportOptions = struct('Color','rg');

numEigenworms = 4; % number of projections to plot

inClusterNeighbourNum = 3;
minNeighbrDist = 1500;
strains = {'npr1','N2'};
wormnums = {'40'}%,'HD'};
analysisTypes = {'loneWorms','inCluster'}%,'smallCluster','leaveCluster'};
pixelsize = 100/19.5; % 100 microns are 19.5 pixels
maxBlobSize = 2.5e5;
minSkelLength = 850;
maxSkelLength = 1500;
postExitDuration = 10;
phase = 'joining';
binwidth = 0.5;

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        figure
        % load information for phase restriction
        [phaseFrames,filenames,~] = xlsread(['datalists/' strain '_' wormnum '_r_list.xlsx'],1,'A1:E15','basic');
        %% load data
        numFiles = length(filenames);
        for analysisType = analysisTypes
            eigenProjections1 = cell(numFiles,1);
            eigenProjections2 = cell(numFiles,1);
            eigenProjections3 = cell(numFiles,1);
            eigenProjections4 = cell(numFiles,1);
            for fileCtr = 1:numFiles
                filename = filenames{fileCtr};
                if exist(strrep(filename,'skeletons','features'),'file')&&exist(filename,'file')
                    features = h5read(strrep(filename,'skeletons','features'),'/features_timeseries');
                    trajData = h5read(filename,'/trajectories_data');
                    blobFeats = h5read(filename,'/blob_features');
                    skelData = h5read(filename,'/skeleton');
                    frameRate = h5readatt(filename,'/plate_worms','expected_fps');
                    %% filter data
                    % filter by blob size and intensity - for non-manually joined trajectories
                    if contains(filename,'55')||contains(filename,'54')
                        intensityThreshold_r = 80;
                    else
                        intensityThreshold_r = 40;
                    end
                    trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                        intensityThreshold_r,maxBlobSize,false,[]);
                    % filter by skeleton length
                    trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)...
                        &filterSkelLength(skelData,pixelsize,minSkelLength,maxSkelLength,false,[]);
                    % apply phase restriction
                    [firstFrame, lastFrame] = getPhaseRestrictionFrames(phaseFrames,phase,fileCtr);
                    phaseFrameLogInd = trajData.frame_number < lastFrame & trajData.frame_number > firstFrame;
                    trajData.filtered = trajData.filtered&phaseFrameLogInd;
                    % filter for in-cluster etc
                    [leaveClusterLogInd, loneWormLogInd, inClusterLogInd,smallClusterLogInd] =...
                        findWormCategory(filename,inClusterNeighbourNum,minNeighbrDist,postExitDuration);
                    switch analysisType{1}
                        case 'loneWorms'
                            features.filtered = ismember(features.skeleton_id+1,find(loneWormLogInd&trajData.filtered));
                        case 'inCluster'
                            features.filtered = ismember(features.skeleton_id+1,find(inClusterLogInd&trajData.filtered));
                        case 'smallCluster'
                            features.filtered = ismember(features.skeleton_id+1,find(smallClusterLogInd&trajData.filtered));
                        case 'leaveCluster'
                            features.filtered = ismember(features.skeleton_id+1,find(leaveClusterLogInd&trajData.filtered));
                    end
                    %% get eigenworm projections
                    if any(features.filtered)
                        eigenProjections1{fileCtr} = features.eigen_projection_1(features.filtered);
                        eigenProjections2{fileCtr} = features.eigen_projection_2(features.filtered);
                        eigenProjections3{fileCtr} = features.eigen_projection_3(features.filtered);
                        eigenProjections4{fileCtr} = features.eigen_projection_4(features.filtered);
                    else
                        warning(['All ' analysisType{1} 'worms filtered out for ' filename ])
                    end
                else
                    warning(['Not all necessary tracking results present for ' filename])
                end
            end
            % pool data from multiple recordings
            eigenProjections1 = vertcat(eigenProjections1{:});
            eigenProjections2 = vertcat(eigenProjections2{:});
            eigenProjections3 = vertcat(eigenProjections3{:});
            eigenProjections4 = vertcat(eigenProjections4{:});
            %% plot data
            subplot(2,2,1)
            if strcmp(analysisType{1},'loneWorms')
                hold on, title(['eigenworm 1'],'FontWeight','normal')
                ylabel('P'),
                yticks(get(gca,'ylim')), xlim([-8 8]),
                xlabel('a_1'), box on
            end
            histogram(eigenProjections1,'Normalization','Probability','DisplayStyle','stairs','BinWidth',binwidth)
            
            subplot(2,2,2)
            if strcmp(analysisType{1},'loneWorms')
                hold on, title(['eigenworm 2'],'FontWeight','normal')
                ylabel('P'),
                yticks(get(gca,'ylim')), xlim([-8 8]),
                xlabel('a_2'), box on
            end
            histogram(eigenProjections2,'Normalization','Probability','DisplayStyle','stairs','BinWidth',binwidth)
            
            subplot(2,2,3)
            if strcmp(analysisType{1},'loneWorms')
                hold on, title(['eigenworm 3'],'FontWeight','normal')
                ylabel('P'),
                yticks(get(gca,'ylim')), xlim([-8 8]),
                xlabel('a_3'), box on
            end
            histogram(eigenProjections3,'Normalization','Probability','DisplayStyle','stairs','BinWidth',binwidth)
            
            subplot(2,2,4)
            if strcmp(analysisType{1},'loneWorms')
                hold on, title(['eigenworm 4'],'FontWeight','normal')
                ylabel('P'),
                yticks(get(gca,'ylim')), xlim([-8 8]) ,
                xlabel('a_4'), box on
            end
            histogram(eigenProjections4,'Normalization','Probability','DisplayStyle','stairs','BinWidth',binwidth)
        end
        %% format and export figure
    legend({'lone worm','in cluster'})
    figFileName = ['figures/eigenwormProjectionDistribution_' strain '_' wormnum '.eps'];
    exportfig(gcf,figFileName,'Color','rgb')
    system(['epstopdf ' figFileName]);
    system(['rm ' figFileName]);
    end
end
%tilefigs()