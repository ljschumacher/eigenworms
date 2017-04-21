clear
% issues/todo:
% - parfor does not properly work for case wormnum 1W

% figure export options
exportOptions = struct('Color','rgb');

% frames to use for calculation of eigenworms, for each combination of
% strain and worm number
numSamples = 10000;

% load master eigenworms for projections
load('singleWorm/masterEigenWorms_N2.mat','eigenWorms');
masterWorms = eigenWorms;

showPlots = true;
plotDiagnostics = true;
nEigenworms = 6;

pixelsize = 100/19.5; % 100 microns are 19.5 pixels
neighbourCutOff = 500; % distance in microns to consider a neighbour close
minNumNeighbours = 3;
strains = {'npr1','N2'};
wormnums = {'40','HD','1W'};
maxBlobSize = 2.5e5;
maxBlobSize_g = 1e4;
minSkelLength = 850;
maxSkelLength = 1500;
intensityThresholds_g = [60, 40, NaN];

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        close all
        %% load data
        % not all results may be present, so check how many
        filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames);
        if ~strcmp(wormnum,'1W')
            filenames_g = importdata(['datalists/' strains{strainCtr} '_' wormnum '_g_list.txt']);
        else
            filenames_g = {};
        end
        % allocate variables
        skeleta = cell(numFiles,1);
        wormIDs = cell(numFiles,1);
        frameIDs = cell(numFiles,1);
        for fileCtr=1:numFiles % can be parfor
            filename = filenames{fileCtr};
            if ~strcmp(wormnum,'1W'), filename_g = filenames_g{fileCtr}; end
            if exist(filename,'file')&&(strcmp(wormnum,'1W')||exist(filename_g,'file'))
                frameRate = h5readatt(filename,'/plate_worms','expected_fps');
                blobFeats = h5read(filename,'/blob_features');
                trajData = h5read(filename,'/trajectories_data');
                skelData = h5read(filename,'/skeleton');
                if ~strcmp(wormnum,'1W')
                    trajData_g = h5read(filename_g,'/trajectories_data');
                    blobFeats_g = h5read(filename_g,'/blob_features');
                end
                %% filter data
                % filter by blob size and intensity
                if contains(filename,'55')||contains(filename,'54')
                    intensityThreshold = 80;
                else
                    intensityThreshold = 40;
                end
                trajData.filtered = filterIntensityAndSize(blobFeats,pixelsize,...
                    intensityThreshold,maxBlobSize,plotDiagnostics,...
                    [strain ' ' wormnum ' ' strrep(filename(end-31:end-5),'/','')]);
                % filter by skeleton length
                trajData.filtered = trajData.filtered&logical(trajData.is_good_skel)...
                    &filterSkelLength(skelData,pixelsize,minSkelLength,maxSkelLength,plotDiagnostics,...
                    [strain '_' wormnum '_' strrep(filename(end-31:end-5),'/','')]);
                % load skeleton data
                if strcmp(wormnum,'1W')
                    skeleta{fileCtr} = skelData(:,:,trajData.filtered);
                else % if it is multiworm data, we need to filter for worms in clusters
                    framesAnalyzed = unique(trajData.frame_number(trajData.filtered));
                    % filter green channel by blob size and intensity
                    trajData_g.filtered = filterIntensityAndSize(blobFeats_g,pixelsize,...
                    intensityThresholds_g(numCtr),maxBlobSize_g);
                    % filter for in-cluster
                    num_close_neighbours_rg = h5read(filename,'/num_close_neighbours_rg');
                    trajData.filtered = num_close_neighbours_rg>=minNumNeighbours;
                    skeleta{fileCtr} = skelData(:,:,trajData.filtered);
                end
            else
                warning(['Not all necessary tracking results present for ' filename ])
            end
        end
        % pool data from multiple recordings
        skeleta = cat(3,skeleta{:});
        % randomly pick numFrames from the data, to not oversample
        if numSamples<=size(skeleta,3)
            [skeleta, framesSampled] = datasample(skeleta,numSamples,3,'Replace',false);
        else
            warning('not enough skeleta to sample form')
        end
        % create angle arrays from skeleta
        [angleArray, ~] = makeAngleArrayV(squeeze(skeleta(1,:,:))',squeeze(skeleta(2,:,:))');
        
        %% find eigenworms
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
        save(['results/eigenData_' strain '_' wormnum '_bodywall.mat'],'eigenWorms',...
            'eigenVals','eigenProjections','varExplained','masterProjections',...
            'masterVarExplained','numFiles','numSamples')
        if showPlots
            % save plots
            figName = [strain '_' wormnum ];
            figPrefix = {'var','eig','cov','eigenValueDistribution'};
            for figCtr = 1:4
                set(figure(figCtr),'name',[figPrefix{figCtr} '_' figName ...
                    '_' num2str(numFiles) 'datasets_' num2str(size(angleArray,1),2) 'frames'])
                figFileName = ['figures/diagnostics/' figPrefix{figCtr} '_' figName '.eps'];
                exportfig(figure(figCtr),figFileName,exportOptions)
                system(['epstopdf ' figFileName]);
                system(['rm ' figFileName]);
            end
        end
    end
end