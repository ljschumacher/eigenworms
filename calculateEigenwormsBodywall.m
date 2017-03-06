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
minNumNeighbours = 2;
strains = {'npr1','N2'};
wormnums = {'1W','40','HD'};
maxBlobSize = 2.5e5;
maxBlobSize_g = 1e4;
minSkelLength = 850;
intensityThresholds_g = [0, 60, 40];

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        close all
        %% load data
        % not all results may be present, so check how many
        filenames = importdata([strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames);
        if ~strcmp(wormnum,'1W')
            filenames_g = importdata([strains{strainCtr} '_' wormnum '_g_list.txt']);
        end
        % allocate variables
        skeleta = cell(numFiles,1);
        wormIDs = cell(numFiles,1);
        frameIDs = cell(numFiles,1);
        for fileCtr=1:numFiles % can be parfor
            filename = filenames{fileCtr};
            frameRate = h5readatt(filename,'/plate_worms','expected_fps');
            blobFeats = h5read(filename,'/blob_features');
            trajData = h5read(filename,'/trajectories_data');
            skelData = h5read(filename,'/skeleton');
            if ~strcmp(wormnum,'1W')
                filename_g = filenames_g{fileCtr};
                trajData_g = h5read(filename_g,'/trajectories_data');
                blobFeats_g = h5read(filename_g,'/blob_features');
            end
            %% filter data
            % filter by blob size and intensity
            if contains(filename,'55')||contains(filename,'54')
                intensityThreshold = 70;
            else
                intensityThreshold = 35;
            end
            trajData.filtered = (blobFeats.area*pixelsize^2<=maxBlobSize)&...
                (blobFeats.intensity_mean>=intensityThreshold)&...
                logical(trajData.is_good_skel);
            % filter by skeleton length
            if plotDiagnostics
                plotSkelLengthDist(skelData(:,:,trajData.filtered),pixelsize,minSkelLength,...
                    [strain '_' wormnum '_' strrep(filename(end-31:end),'/','')]);
            end
            trajData.filtered(...
                sum(sqrt(sum((diff(skelData,1,2)*pixelsize).^2)))<minSkelLength)...
                = false;
            % load skeleton data
            if strcmp(wormnum,'1W')
                skeleta{fileCtr} = skelData(:,:,trajData.filtered);
            else % if it is multiworm data, we need to filter for worms in clusters
                framesAnalyzed = unique(trajData.frame_number(trajData.filtered));
                numFrames = numel(framesAnalyzed);
                % filter green channel by blob size and intensity
                trajData_g.filtered = (blobFeats_g.area*pixelsize^2<=maxBlobSize_g)&...
                    (blobFeats_g.intensity_mean>=intensityThresholds_g(numCtr));
                % calculate red-green neighbour distances to filter for in-cluster
                for frameCtr = 1:numFrames
                    frame = framesAnalyzed(frameCtr);
                    [x, y] = getWormPositions(trajData, frame);
                    [x_g, y_g] = getWormPositions(trajData_g, frame);
                    if numel(x_g)>=1&&numel(x)>=1 % need at least two worms in frame to calculate distances
                        redToGreenDistances = pdist2([x y],[x_g y_g]).*pixelsize; % distance of every red worm to every green
                        numNeighbours = sum(redToGreenDistances<neighbourCutOff,2);
                    elseif numel(x)>=1
                        numNeighbours = zeros(size(x));
                    end
                    % exclude worms with fewer than some number of neighbours
                    trajData.filtered(...
                        (trajData.frame_number==frame)&trajData.filtered) = ...
                        numNeighbours>=minNumNeighbours; % keeps only those worms more than minNumNeighbours neighbours
                end
                skeleta{fileCtr} = skelData(:,:,trajData.filtered);
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