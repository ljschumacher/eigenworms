clear

% figure export options
exportOptions = struct('Color','rgb');

% frames to use for calculation of eigenworms, for each combination of
% strain and worm number
numFrames = 10000;

% load master eigenworms for projections
load('singleWorm/masterEigenWorms_N2.mat','eigenWorms');
masterWorms = eigenWorms;

showPlots = 1;
nEigenworms = 6;

pixelsize = 100/19.5; % 100 microns are 19.5 pixels

strains = {'npr1','N2'};
wormnums = {'40','HD','1W'};
maxBlobSize = 2.5e5;

for numCtr = 1:length(wormnums)
    wormnum = wormnums{numCtr};
    for strainCtr = 1:length(strains)
        strain = strains{strainCtr};
        close all
        %% load data
        % not all results may be present, so check how many
        filenames = importdata([strains{strainCtr} '_' wormnum '_r_list.txt']);
        numFiles = length(filenames);
        % allocate variables
        skeleta = cell(numFiles,1);
        wormIDs = cell(numFiles,1);
        frameIDs = cell(numFiles,1);
        skeletaSizes = NaN(numFiles,1);
        for fileCtr=1:numFiles
            filename = filenames{fileCtr};
            frameRate = h5readatt(filename,'/plate_worms','expected_fps');
            %             blobFeats = h5read(filename,'/blob_features');
            trajData = h5read(filename,'/trajectories_data');
            skelData = h5read(filename,'/skeleton');
            %% filter data
            if strcmp(wormnum,'1W')
                skeleta{fileCtr} = skelData(:,:,logical(trajData.is_good_skel));
            else
                % if it is multiworm data, we need to filter for isolated worms
                loneWorms = horzcat(mindist{:})>=2500;
                skeleta{fileCtr} = skelData(:,:,loneWorms&logical(trajData.is_good_skel));
            end
        end
        % pool data from multiple recordings
        skeleta = cat(3,skeleta{:});
        % randomly pick numFrames from the data, to not oversample
        [skeleta, framesSampled] = datasample(skeleta,numFrames,3,'Replace',false);
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
            'masterVarExplained','numFiles','numFrames')
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