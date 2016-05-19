close all
clear

% figure export options
exportOptions = struct('Color','rgb');

% select data set by strain - N2, HW, NP
for strain = {'N2', 'HW', 'NP'}
    S = strain{:};
    % select data set by number of worms - 1, 5, 15, 25, 40
    for N = [40 25 15 5 1]
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
            
            for dsCtr=1:nDatasets
                filename = datasets{dsCtr};
                % load trajectory data
                trajectoryData = h5read(filename,'/trajectories_data');
                % plot worm numbers
                frameNums = unique(trajectoryData.frame_number)';
                histogram(trajectoryData.frame_number,frameNums,'DisplayStyle','stairs')
                hold on
                % filter data
                % select frames that have skeleton
                hasSkel = trajectoryData.has_skeleton==1;
                % select frames with worms that occurr more than a certain number
                framesPerWorm = histcounts(trajectoryData.worm_index_joined,max(trajectoryData.worm_index_joined));
                frequentWorms = find(framesPerWorm>=25*30);
                framesFilter = ismember(trajectoryData.worm_index_joined,frequentWorms);
                % select frames with a certain area
                areaFilter = filterArea(trajectoryData,25,1500,50,50,hasSkel&framesFilter);
                % detect dust from manually labelled data
                dustIdcs = filterDust(trajectoryData,0.1,5);
                % select worms with at least a certain speed
                speedFilter = filterSpeed(trajectoryData,0.1,1,hasSkel&framesFilter&areaFilter&~dustIdcs);
                % combine filters
                combiFilter = hasSkel&areaFilter&framesFilter&speedFilter&~dustIdcs;
                % plot post-filtered worm numbers
                filteredFrameNums = unique(trajectoryData.frame_number(combiFilter))';
                histogram(trajectoryData.frame_number(combiFilter),filteredFrameNums,'DisplayStyle','stairs')
                xlabel('frame number')
                ylabel('object count')
                legend('raw','filtered')
                % save plot
                figName = [S '_' num2str(N) 'worms' '_dataset_' filename(end-42:end-15)];
                figSuffix = {'var','eig','cov'};
                    set(gcf,'name',figName)
                    figFileName = ['figures/diagnostics' figName '_wormNums.eps'];
                    exportfig(gcf,figFileName,exportOptions)
                    system(['epstopdf ' figFileName]);
                    system(['rm ' figFileName]);
                hold off
            end
        else
            display(['No datasets found for strain=' S ', worms=' num2str(N) ])
        end
    end
end
% histogram(areas)
tilefigs%([5 6])