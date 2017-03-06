% combine projected amplitudes in space of first two eigenworms

% issues/to-do:
% -
close all
clear

% figure export options
exportOptions = struct('Color','rgb');

nComponents = 4;

strains = {'N2', 'npr1'};
nStrains = length(strains);
wormnums = {'HD','40','1W'};
for numCtr = 1:length(wormnums)
    N = wormnums{numCtr};
    eigProjectionFig = figure;
    % loop through different strains
    for strainCtr = 1:nStrains
        S = strains{strainCtr};
        % load eigenworm data
        file = rdir(['results/eigenData_' S '_' N '_bodywall.mat']);
        if ~isempty(file)
            % load eigenworm analysis result
            load(file.name,'masterProjections');
%             % change the order of first and third masterworm
%             masterProjections = masterProjections(:,[3 2 1 4:end]);
            % normalise to unit variance
% %             masterProjections = zscore(masterProjections);
            % plot projected amplitudes
            for cmpCtr = 1:nComponents
                subplot(ceil(sqrt(nComponents)),floor(sqrt(nComponents)),cmpCtr)
                histogram(masterProjections(:,cmpCtr),...
                    'Normalization','Probability','DisplayStyle','stairs')
                if strainCtr == 1
                    hold on
                    ax = gca;
                    ax.YLabel.String = 'P';
                    ax.XLabel.String =  ['a_' num2str(cmpCtr)];
                    ax.XLim = [-10 10];
                end
            end
            clear masterProjections
        else
            display(['No data for strain=' S ', worms=' num2str(N)])
        end
    end
    % annotate and save figure
    legend(strains,'Location','East')
    set(eigProjectionFig, 'name', ['projected amplitudes ' N ' worm data'])
    figName = ['figures/projections_' N 'worms_strainComparison.eps'];
    exportfig(eigProjectionFig,figName,exportOptions)
    system(['epstopdf ' figName]);
    system(['rm ' figName]);
    %             close(eigProjectionFig)
end
tilefigs([2 3])