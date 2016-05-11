% combine projected amplitudes in space of first two eigenworms

% issues/to-do:
% - should we normalise the variance of each projection to 1, as in the
% original eigenworm paper?
close all
clear

% figure export options
exportOptions = struct('Color','rgb','Renderer','painters');

nComponents = 4;

strains = {'N2', 'HW', 'NP'};
nStrains = length(strains);
for N = [1 5 15 25 40]
    % loop through different strains
    for strainCtr = 1:nStrains
        S = strains{strainCtr};
        % load eigenworm data
        file = rdir(['results/' S '_' num2str(N) 'worms_eigenData.mat']);
        if ~isempty(file)
            eigProjectionFig = figure;
            % load eigenworm analysis result
            load(file.name)
            % normalise to unit variance
            eigenProjections = zscore(eigenProjections);
            % plot projected amplitudes
            [~, subAx, bigAx, ~, diagAx] = hplotmatrix(eigenProjections(:,1:nComponents));
            for cmpCtr = 1:nComponents
                subAx(cmpCtr).YLabel.String = ['a_' num2str(cmpCtr)];
                subAx(cmpCtr*nComponents).XLabel.String = ['a_' num2str(cmpCtr)];
                diagAx(cmpCtr).XLim = [-2 2];
                for plotCtr = 1:nComponents
                    ax = subAx(cmpCtr,plotCtr);
                    ax.XLim = [-2 2]; ax.YLim = [-2 2];
                end
            end
            bigAx.Title.String = ['strain ' S ', N = ' num2str(N) ' worms, '...
                num2str(size(eigenProjections,1)/25/3600,2) ' worm-hours'];
            % annotate and save figure
            set(eigProjectionFig, 'name', ['projected amplitudes for N=' num2str(N) ' worms'])
            figName = ['figures/' S '_' num2str(N) 'worms_projections.eps'];
            exportfig(eigProjectionFig,figName,exportOptions)
            system(['epstopdf ' figName]);
            system(['rm ' figName]);
            %             close(eigProjectionFig)
            clear eigenProjections
        else
            display(['No data for strain=' S ', worms=' num2str(N)])
        end
    end
end
tilefigs([3 5])