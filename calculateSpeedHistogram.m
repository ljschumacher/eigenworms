close all
clear

areas = [];

dBin = 0.05;
bins = 0:dBin:2.5;
strains = {'N2', 'HW', 'NP'};
% select data set by strain - N2, HW, NP
for strainCtr = 1:length(strains)
    S = strains{strainCtr};
    % select data set by number of worms - 1, 5, 15, 25, 40
    for N = [1 5 15 25 40]
        % %         close all
        
        % load file name descriptor - taken from Camille's Recording_LOG.xls
        filelist = load(['recordingsLog/strain' S 'worms' num2str(N) '.mat']);
        filenames = filelist.filenames;
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
            
            parfor dsCtr=1:nDatasets
                filename = datasets{dsCtr};
                % load trajectory data
                trajectoryData = h5read(filename,'/trajectories_data');
                hasSkel = trajectoryData.has_skeleton==1;
                framesPerWorm = histcounts(trajectoryData.worm_index_joined,max(trajectoryData.worm_index_joined));
                frequentWorms = find(framesPerWorm>=25*30);
                frequentFilter = ismember(trajectoryData.worm_index_joined,frequentWorms);
                % select frames with a certain area
                areaFilter = filterArea(trajectoryData,20,50,hasSkel&frequentFilter,0,filename);
                % find individual worm trajectories
                wormIDs = unique(trajectoryData.worm_index_joined(hasSkel&frequentFilter&areaFilter))';
                % go through worms and calculate speeds
                speedHists = NaN(length(wormIDs),length(bins)-1);
                for wormCtr=1:length(wormIDs)
                    wormID = wormIDs(wormCtr);
                    wormIdcs = trajectoryData.worm_index_joined==wormID;
                    wormDx = diff(trajectoryData.coord_x(wormIdcs));
                    wormDy = diff(trajectoryData.coord_y(wormIdcs));
                    wormDf = diff(trajectoryData.frame_number(wormIdcs));
                    wormSpeed = sqrt(wormDx.^2 + wormDy.^2)./single(wormDf);
                    speedHists(wormCtr,:) = histcounts(wormSpeed,bins,'Normalization','pdf');
                end
                figure, plot(bins(1:end-1)+dBin/2,speedHists')
                title([S ' N=' num2str(N)])
                set(gcf,'name',filename)
            end
        else
            display(['No datasets found for strain=' S ', worms=' num2str(N) ])
        end
    end
end
% histogram(areas)
tilefigs%([5 6])