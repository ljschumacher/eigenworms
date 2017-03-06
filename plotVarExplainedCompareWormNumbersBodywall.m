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
plotColors = lines(length(wormnums));
for strainCtr = 1:nStrains
    S = strains{strainCtr};
    varFig = figure;
    plotHandles = NaN(size(wormnums));
    for numCtr = 1:length(wormnums)
        N = wormnums{numCtr};
        % load eigenworm data
        file = rdir(['results/eigenData_' S '_' N '_bodywall.mat']);
        if ~isempty(file)
            % load eigenworm analysis result
            load(file.name,'masterVarExplained','varExplained');
%             % change the order of first and third masterworm to match ours
%             masterVarExplained = masterVarExplained([3 2 1 4:end]);
            plotHandles(numCtr) =...
                plot(0:6,cumsum([0; varExplained]),'-o','Color',plotColors(numCtr,:));            
            if numCtr==1
                hold on
                ax = gca;
                ax.YLabel.String = '% variance explained';
                ax.XLabel.String = 'eigenworm';
            end
            plot(0:6,cumsum([0; masterVarExplained]),'--x','Color',plotColors(numCtr,:))
        else
            display(['No data for strain=' S ', worms=' num2str(N)])
        end
    end
    % annotate and save figure
    legend(plotHandles, wormnums,'Location','SouthEast')
    title([S ' variance explained by own (-o) and reference (--x) eigenworms'],'FontWeight','normal')
    figName = ['figures/varExplained_' S '_CompareWormNumbers.eps'];
    exportfig(varFig,figName,exportOptions)
    system(['epstopdf ' figName]);
    system(['rm ' figName]);
    %             close(eigProjectionFig)
end
tilefigs([3 5])