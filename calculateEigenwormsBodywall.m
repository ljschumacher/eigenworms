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
plotDiagnostics = false;
nEigenworms = 6;

pixelsize = 100/19.5; % 100 microns are 19.5 pixels
neighbrCutOff = 500; % distance in microns to consider a neighbr close
minNumNeighbrs = 3;
strains = {'npr1','N2'};
wormnums = {'40','HD','1W'};
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
        filenames = importdata(['datalists/' strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames);
        % allocate variables
        skeleta = cell(numFiles,1);
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
                    % filter for in-cluster
                    num_close_neighbrs_rg = h5read(filename,'/num_close_neighbrs_rg')';
                    trajData.filtered = trajData.filtered&num_close_neighbrs_rg>=minNumNeighbrs;
                    skeleta{fileCtr} = skelData(:,:,trajData.filtered);
                    assert(~any(isnan(skeleta{fileCtr}(:))))
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