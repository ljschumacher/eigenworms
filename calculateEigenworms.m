clear

% select data set by strain - N2, HW, NP
for strain = {'N2', 'HW', 'NP'}
    S = strain{:};
    % select data set by number of worms - 1, 5, 15, 25, 40
    for N = [1 5 15 25 40]
        close all
        % select number of (consecutive) frames to load from each dataset
        framesPerRead = 100000;
        
        % load file name descriptor - taken from Camille's Recording_LOG.xls
        load(['recordingsLog/strain' S 'worms' num2str(N) '.mat']);
        nFiles = length(filenames);
        % not all results may be present, so check how many
        dsCtr = 1;
        datasets = {};
        for fileCtr=1:nFiles
            % find full path to folder
            file = rdir(['/data1/linus/Recordings/Results/*/' ...
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
            skeletaSizes = NaN(nDatasets,1);
            
            for dsCtr=1:nDatasets
                filename = datasets{dsCtr};
                % load skeleton data in chunks of nFrames
                skeleta{dsCtr} = filterData(filename,framesPerRead,1);
                skeletaSizes(dsCtr) = size(skeleta{dsCtr},3);
            end
            
            nSkeleta = sum(skeletaSizes);
            angleArray = NaN(nSkeleta,48);
            for dsCtr=1:nDatasets
                % create angle arrays from skeleta
                [angleArray(sum(skeletaSizes(1:(dsCtr - 1))) + (1:skeletaSizes(dsCtr)),:), ~] = ...
                    makeAngleArrayV(squeeze(skeleta{dsCtr}(1,:,:))',squeeze(skeleta{dsCtr}(2,:,:))');
            end
            
            % find eigenworms
            showPlots = 1;
            nEigenworms = 6;
            [eigenWorms, eigenVals] = findEigenWorms(angleArray, nEigenworms, showPlots);
            eigenProjections = projectOnEigenWorms(eigenWorms, angleArray, nEigenworms);
            
            % save eigenWorms, eigenVals and first few projections
            nFrames = num2str(size(angleArray,1),2);
            save(['results/' S '_' num2str(N) 'worms_eigenData.mat'],'eigenWorms',...
                'eigenVals','eigenProjections','nDatasets','nFrames')
            if showPlots
                % save plots
                figName = [S '_' num2str(N) 'worms_' num2str(nDatasets) 'datasets_' ...
                    num2str(nFrames,2) 'frames'];
                figSuffix = {'var','eig','cov'};
                for figCtr = 1:3
                    set(figure(figCtr),'name',[figName '_' figSuffix{figCtr}])
                    savefig(figure(figCtr),['figures/' S '_' num2str(N) ...
                        'worms_' figSuffix{figCtr} '.fig'],'compact')
                end
            end
        else
            display(['No datasets found for strain=' S ', worms=' num2str(N) ])
        end
    end
end