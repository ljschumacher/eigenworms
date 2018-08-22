% combine projected amplitudes in space of first two eigenworms

% issues/to-do:
% - 
close all
clear

% figure export options
exportOptions = struct('Color','rgb','Renderer','painters');

nComponents = 4;

strains = {'N2', 'npr1'};
nStrains = length(strains);
wormnums = {'40'}%{'HD','40','1W'};
for numCtr = 1:length(wormnums)
    N = wormnums{numCtr};
    % loop through different strains
    for strainCtr = 1:nStrains
        S = strains{strainCtr};
        % load eigenworm data
        file = rdir(['results/eigenData_' S '_' N '_bodywall.mat']);
        if ~isempty(file)
            eigProjectionFig = figure;
            % load eigenworm analysis result
            load(file.name,'masterProjections');
            % change the order of first and third masterworm to match ours
            masterProjections = masterProjections(:,[3 2 1 4:end]);
            % normalise to unit variance
            masterProjections = zscore(masterProjections);
            % plot projected amplitudes
            [~, subAx, bigAx, ~, diagAx] = hplotmatrix(masterProjections(:,1:nComponents));
            for cmpCtr = 1:nComponents
                subAx(cmpCtr).YLabel.String = ['a_' num2str(cmpCtr)];
                subAx(cmpCtr*nComponents).XLabel.String = ['a_' num2str(cmpCtr)];
                diagAx(cmpCtr).XLim = [-2 2];
                for plotCtr = 1:nComponents
                    ax = subAx(cmpCtr,plotCtr);
                    ax.XLim = [-2 2]; ax.YLim = [-2 2];
                end
            end
            bigAx.Title.String = ['strain ' S ', N = ' N ' worms, '...
                num2str(size(masterProjections,1)) ' worm-frames'];
            % annotate and save figure
            set(eigProjectionFig, 'name', ['projected amplitudes for N=' N ' worms'])
            figName = ['figures/projections_' S '_' N 'worms.eps'];
            exportfig(eigProjectionFig,figName,exportOptions)
            system(['epstopdf ' figName]);
            system(['rm ' figName]);
            %             close(eigProjectionFig)
            clear masterProjections
        else
            display(['No data for strain=' S ', worms=' num2str(N)])
        end
    end
end
tilefigs([3 5])