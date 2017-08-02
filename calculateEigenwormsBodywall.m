clear
% issues/todo:
% - parfor does not properly work for case wormnum 1W

% specify whether to phase-restrict
phase = 'fullMovie'; % 'fullMovie' or 'stationary'

% figure export options
exportOptions = struct('Color','rgb');

% frames to use for calculation of eigenworms, for each combination of
% strain and worm number
numSamples = 10000;

% load master eigenworms for projections
load('singleWorm/masterEigenWorms_N2.mat','eigenWorms');
masterWorms = eigenWorms;

showPlots = true;
plotDiagnostics = false;
nEigenworms = 6;

pixelsize = 100/19.5; % 100 microns are 19.5 pixels
neighbrCutOff = 500; % distance in microns to consider a neighbr close
minNumNeighbrs = 3;
minNeighbrDist = 1500;
strains = {'npr1','N2'};
wormnums = {'1W','40','HD'};
if strcmp(phase, 'stationary')
    wormnums = {'40'};
end

maxBlobSize = 2.5e5;
minSkelLength = 850;
maxSkelLength = 1500;

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        close all
        %% load data
        % not all results may be present, so check how many
        if strcmp(phase, 'stationary')
            [lastFrames,filenames,~] = xlsread(['../trackingAnalysis/datalists/' strains{strainCtr} '_' wormnum '_r_list.xlsx'],1,'A1:B15','basic');
        else
            filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        end
        numFiles = length(filenames);
        % allocate variables
        skeletaLoneWorms = cell(numFiles,1);
        if ~strcmp(wormnum,'1W')
            skeletaInCluster = cell(numFiles,1);
            skeletaSmallCluster = cell(numFiles,1);
            skeletaLeaveCluster = cell(numFiles,1);
        end
        wormIDs = cell(numFiles,1);
        frameIDs = cell(numFiles,1);
        for fileCtr=1:numFiles % can be parfor
            filename = filenames{fileCtr};
            if exist(filename,'file')
                frameRate = h5readatt(filename,'/plate_worms','expected_fps');
                blobFeats = h5read(filename,'/blob_features');
                trajData = h5read(filename,'/trajectories_data');
                skelData = h5read(filename,'/skeleton');
                %% filter data
                % filter by blob size and intensity
                if contains(filename,'55')||contains(filename,'54')
                    intensityThreshold_r = 80;
                else
                    intensityThreshold_r = 40;
                end
                trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                    intensityThreshold_r,maxBlobSize,plotDiagnostics,...
                    [strain ' ' wormnum ' ' strrep(filename(end-31:end-5),'/','')]);
                % filter by skeleton length
                trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)...
                    &filterSkelLength(skelData,pixelsize,minSkelLength,maxSkelLength,plotDiagnostics,...
                    [strain '_' wormnum '_' strrep(filename(end-31:end-5),'/','')]);
                % load skeleton data
                if strcmp(wormnum,'1W')
                    skeletaLoneWorms{fileCtr} = skelData(:,:,trajData.filtered);
                else % if it is multiworm data, we need to filter for worms in clusters
                    % filter for in-cluster etc
                    min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
                    num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
                    neighbr_dist = h5read(filename,'/neighbr_distances');
                    loneWorms = min_neighbr_dist>=minNeighbrDist;
                    inCluster = num_close_neighbrs>=minNumNeighbrs;
                    smallCluster = (num_close_neighbrs==2 & neighbr_dist(:,3)>=minNeighbrDist)...
                        |(num_close_neighbrs==3 & neighbr_dist(:,4)>=minNeighbrDist)...
                        |(num_close_neighbrs==4 & neighbr_dist(:,5)>=minNeighbrDist);
                    leaveCluster = [false; inCluster(1:end-1)&~inCluster(2:end)];
                    leaveClusterExtended = unique(find(leaveCluster) + [0:5*double(frameRate)]);% add 5 more second after each leaving event
                    leaveClusterExtended = leaveClusterExtended(leaveClusterExtended<numel(leaveCluster)); % exclude frames beyond highest frame number
                    leaveCluster(leaveClusterExtended) = true;
                    leaveCluster(inCluster) = false; % exclude times when worm moved back
                    leaveCluster(loneWorms)=false; % exclude worms that have become lone worm
                    % restrict analysis to stationary phase only
                    if strcmp(phase,'stationary')
                        firstFrame = double(round(max(trajData.frame_number)/10)); % define starting frame as 10% into the movie
                        lastFrame = lastFrames(fileCtr);
                        phaseFrameLogInd = trajData.frame_number < lastFrame & trajData.frame_number > firstFrame;
                        trajData.filtered(~phaseFrameLogInd) = false;
                    end
                    % load skeletal data and check it's not NaN
                    skeletaLoneWorms{fileCtr} = skelData(:,:,trajData.filtered&loneWorms);
                    skeletaInCluster{fileCtr} = skelData(:,:,trajData.filtered&inCluster);
                    skeletaSmallCluster{fileCtr} = skelData(:,:,trajData.filtered&smallCluster);
                    skeletaLeaveCluster{fileCtr} = skelData(:,:,trajData.filtered&leaveCluster);
                    assert(~any(isnan(skeletaInCluster{fileCtr}(:))))
                    assert(~any(isnan(skeletaSmallCluster{fileCtr}(:))))
                    assert(~any(isnan(skeletaLeaveCluster{fileCtr}(:))))
                end
                assert(~any(isnan(skeletaLoneWorms{fileCtr}(:))))
            else
                warning(['Not all necessary tracking results present for ' filename ])
            end
        end
        % pool data from multiple recordings
        skeletaLoneWorms = cat(3,skeletaLoneWorms{:});
        if ~strcmp(wormnum,'1W')
            skeletaInCluster = cat(3,skeletaInCluster{:});
            skeletaSmallCluster = cat(3,skeletaSmallCluster{:});
            skeletaLeaveCluster = cat(3,skeletaLeaveCluster{:});
        end
        % randomly pick numFrames from the data, to not oversample, and create angle arrays from skeleta
        if numSamples<=size(skeletaLoneWorms,3)
            [skeletaLoneWorms, ~] = datasample(skeletaLoneWorms,numSamples,3,'Replace',false);
        else
            warning(['not enough lone worm skeleta to sample from for ' wormnum ' ' strain])
        end
        [angleArrayLoneWorms, ~] = makeAngleArrayV(squeeze(skeletaLoneWorms(1,:,:))',squeeze(skeletaLoneWorms(2,:,:))');
        clear skeletaLoneWorms % free some memory
        if ~strcmp(wormnum,'1W')
            if numSamples<=size(skeletaInCluster,3)
                [skeletaInCluster, ~] = datasample(skeletaInCluster,numSamples,3,'Replace',false);
            else
                warning(['not enough in cluster worm skeleta to sample from for ' wormnum ' ' strain])
            end
            [angleArrayInCluster, ~] = makeAngleArrayV(squeeze(skeletaInCluster(1,:,:))',squeeze(skeletaInCluster(2,:,:))');
            clear skeletaInCluster % free some memory
            if numSamples<=size(skeletaSmallCluster,3)
                [skeletaSmallCluster, ~] = datasample(skeletaSmallCluster,numSamples,3,'Replace',false);
            else
                warning(['not enough small cluster skeleta to sample from for ' wormnum ' ' strain])
            end
            [angleArraySmallCluster, ~] = makeAngleArrayV(squeeze(skeletaSmallCluster(1,:,:))',squeeze(skeletaSmallCluster(2,:,:))');
            clear skeletaSmallCluster % free some memory
            if numSamples<=size(skeletaLeaveCluster,3)
                [skeletaLeaveCluster, ~] = datasample(skeletaLeaveCluster,numSamples,3,'Replace',false);
            else
                warning(['not enough in cluster worm skeleta to sample from for ' wormnum ' ' strain])
            end
            [angleArrayLeaveCluster, ~] = makeAngleArrayV(squeeze(skeletaLeaveCluster(1,:,:))',squeeze(skeletaLeaveCluster(2,:,:))');
            clear skeletaLeaveCluster % free some memory
        end
        %% find eigenworms
        if ~strcmp(wormnum,'1W')
            analysisTypes = {'loneWorms','inCluster','smallCluster','leaveCluster'};
        else
            analysisTypes = {'loneWorms'};
        end
        for analysisType = analysisTypes
            switch analysisType{1}
                case 'loneWorms'
                    angleArray = angleArrayLoneWorms;
                case 'inCluster'
                    angleArray = angleArrayInCluster;
                case 'smallCluster'
                    angleArray = angleArraySmallCluster;
                case 'leaveCluster'
                    angleArray = angleArrayLeaveCluster;
            end
            [eigenWorms, eigenVals] = findEigenWorms(angleArray, nEigenworms, showPlots);
            % save projections onto reduced dimensions, also for reference
            % components and calc variance explained by reference comps
            eigenProjections = projectOnEigenWormsV(eigenWorms, angleArray, nEigenworms);
            varExplained = eigenVals/sum(var(angleArray));
            masterProjections = projectOnEigenWormsV(masterWorms, angleArray, nEigenworms);
            masterEigVals = diag(masterWorms(1:nEigenworms,:)*cov(angleArray,0,'omitrows')...
                /masterWorms(1:nEigenworms,:));
            masterVarExplained = masterEigVals/sum(var(angleArray));
            % save eigenWorms, eigenVals and first few projections
            save(['results/eigenData_' strain '_' wormnum '_bodywall_' analysisType{1} '_' phase '.mat'],'eigenWorms',...
                'eigenVals','eigenProjections','varExplained','masterProjections',...
                'masterVarExplained','numFiles','numSamples')
            if showPlots
                % save plots
                figName = [strain '_' wormnum '_' phase ];
                figPrefix = {'var','eig','cov','eigenValueDistribution'};
                for figCtr = 1:4
                    set(figure(figCtr),'name',[figPrefix{figCtr} ' ' figName ...
                        ' '  analysisType{1} ' ' num2str(numFiles) ' datasets ' num2str(size(angleArray,1),2) ' frames'])
                    figFileName = ['figures/diagnostics/' figPrefix{figCtr} '_' analysisType{1} '_' figName '.eps'];
                    exportfig(figure(figCtr),figFileName,exportOptions)
                    system(['epstopdf ' figFileName]);
                    system(['rm ' figFileName]);
                end
            end
        end
    end
end