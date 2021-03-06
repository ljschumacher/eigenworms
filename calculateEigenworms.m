clear

% figure export options
exportOptions = struct('Color','rgb');

% frames to use for calculation of eigenworms, for each combination of
% strain and worm number
nFrames = 450000;

% load master eigenworms for projections
load('singleWorm/masterEigenWorms_N2.mat','eigenWorms');
masterWorms = eigenWorms;

showPlots = 1;
nEigenworms = 6;
            
% select data set by strain - N2, HW, NP
for strain = {'N2', 'HW', 'NP'}
    S = strain{:};
    % select data set by number of worms - 1, 5, 15, 25, 40
    for N = [1 5 15 25 40]
        close all
        
        % load file name descriptor - taken from Camille's Recording_LOG.xls
        load(['recordingsLog/strain' S 'worms' num2str(N) '.mat']);
        nFiles = length(filenames);
        % not all results may be present, so check how many
        dsCtr = 0;
        datasets = {};
        for fileCtr=1:nFiles
            % find full path to folder
            file = rdir(['/data1/linus/CamilleData/Results/*/' ...
                filenames{fileCtr}(1:end-5) '_skeletons.hdf5']);
            if ~isempty(file)
                datasets{dsCtr} = file.name;
                dsCtr = dsCtr + 1;
            end
        end
        
        if ~isempty(datasets)
            nDatasets = length(datasets);
            
            % allocate variables
            skeleta = cell(nDatasets,1);
            wormIDs = cell(nDatasets,1);
            frameIDs = cell(nDatasets,1);
            skeletaSizes = NaN(nDatasets,1);
            
            parfor dsCtr=1:nDatasets
                filename = datasets{dsCtr};
                % load skeleton data in chunks of nFrames
                [skeleta{dsCtr}, wormIDs{dsCtr}, frameIDs{dsCtr}] = filterData(filename,1,1);
                skeletaSizes(dsCtr) = size(skeleta{dsCtr},3);
            end
            
            nSkeleta = sum(skeletaSizes);
            angleArray = NaN(nSkeleta,48);
            wormIDarray = int32(zeros(nSkeleta,1,'single'));
            frameIDarray = int32(zeros(nSkeleta,1,'single'));
            for dsCtr=1:nDatasets
                % create angle arrays from skeleta
                [angleArray(sum(skeletaSizes(1:(dsCtr - 1))) + (1:skeletaSizes(dsCtr)),:), ~] = ...
                    makeAngleArrayV(squeeze(skeleta{dsCtr}(1,:,:))',squeeze(skeleta{dsCtr}(2,:,:))');
                % to avoid duplicating worm and frame IDs, pre-append data
                % set number multiplied with next highest order of mag
                wormIDarray(sum(skeletaSizes(1:(dsCtr - 1))) + (1:skeletaSizes(dsCtr)),:) = ...
                    int32(dsCtr*10^ceil(log10(nSkeleta))) + wormIDs{dsCtr};
                frameIDarray(sum(skeletaSizes(1:(dsCtr - 1))) + (1:skeletaSizes(dsCtr)),:) = ...
                    int32(dsCtr*10^ceil(log10(nSkeleta))) + frameIDs{dsCtr};
            end
            
%             display('Calculating shape correlations...')
%             plotShapeCorrelations(angleArray,wormIDarray,frameIDarray,250,...
%                 ['figures/diagnostics/' S '_' num2str(N) 'worms_'],0)
%             masterProjections = projectOnEigenWormsV(masterWorms, angleArray, nEigenworms);
%             plotShapeCorrelations(masterProjections,wormIDarray,frameIDarray,250,...
%                 ['figures/diagnostics/' S '_' num2str(N) 'worms_eigenProj_'],1)

            % randomly pick nFrames from the data, to not oversample
            % from correlated frames
            if nSkeleta >= nFrames
                [angleArray, sampleIDs] = datasample(angleArray,nFrames,'Replace',false);
                wormIDarray = wormIDarray(sampleIDs);
                frameIDarray = frameIDarray(sampleIDs);
                % if not enough frames exist, just take all there are
            end
                        
            % find eigenworms
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
            save(['results/' S '_' num2str(N) 'worms_eigenData.mat'],'eigenWorms',...
                'eigenVals','eigenProjections','varExplained','masterProjections',...
                'masterVarExplained','nDatasets','nFrames','wormIDarray','frameIDarray')
            if showPlots
                % save plots
                figName = [S '_' num2str(N) 'worms'];
                figSuffix = {'var','eig','cov','eigenValueDistribution'};
                for figCtr = 1:4
                    set(figure(figCtr),'name',[figName '_' figSuffix{figCtr} ...
                        '_' num2str(nDatasets) 'datasets_' num2str(size(angleArray,1),2) 'frames'])
                    figFileName = ['figures/diagnostics/' figName '_' figSuffix{figCtr} '.eps'];
                    exportfig(figure(figCtr),figFileName,exportOptions)
                    system(['epstopdf ' figFileName]);
                    system(['rm ' figFileName]);
                end
            end
        else
            display(['No datasets found for strain=' S ', worms=' num2str(N) ])
        end
    end
end