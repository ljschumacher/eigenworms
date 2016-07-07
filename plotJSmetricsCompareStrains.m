% calculate JS metric between distributions of eigencomponents to
% distinguish different strains

close all
clear

% figure export options
exportOptions = struct('Color','rgb','Width',15,'FontMode','scaled','FontSizeMin',1);

nComponents = 4;

binLimits = [-5 5];
binWidth = 0.1;
nBins = length(min(binLimits):binWidth:max(binLimits));
strains = {'N2', 'HW', 'NP'};
textalignments = {'right','left','center'};
lineHandles = NaN(nComponents,1);
nStrains = length(strains);
distances = NaN(nStrains*(nStrains - 1)/2,nComponents);
plotColors = lines(nComponents);
for N = [40]
    distributions =  NaN(nStrains,nComponents,nBins - 1);
    variances = NaN(nStrains,nComponents);
    jsFig = figure;
    hold on
    % loop through different strains
    for strainCtr = 1:nStrains
        S = strains{strainCtr};
        % load eigenworm data
        file = rdir(['results/' S '_' num2str(N) 'worms_eigenData.mat']);
        if ~isempty(file)
            % load eigenworm analysis result
            load(file.name,'masterProjections','masterVarExplained')
            % change the order of first and third masterworm to match ours
            masterProjections = masterProjections(:,[3 2 1 4:nComponents]);
            variances(strainCtr,:) = masterVarExplained([3 2 1 4:nComponents]);
            % normalise to unit variance
            masterProjections = zscore(masterProjections);
            % loop through projected components
            for dimCtr = 1:nComponents
                % calculate 1D histograms
                distributions(strainCtr,dimCtr,:) = histcounts(masterProjections(:,dimCtr),...
                    'BinWidth',binWidth,'BinLimits',binLimits,'Normalization','Probability');
            end
            clear masterProjections
        else
            display(['No data for strain=' S ', worms=' num2str(N)])
        end
    end
    % for each component calculate js-metric between all pairs of strains
    for dimCtr = 1:nComponents
        for i = 1:nStrains
            for j = (i + 1):nStrains
                pairCtr = (i - 1)*(nStrains - i/2) + j - i;
                distances(pairCtr,dimCtr) = ...
                    sqrt(jsdiv(distributions(i,dimCtr,:),distributions(j,dimCtr,:)));
                % plot line at height=distance and length = difference in
                % variance explained
                lineHandles(dimCtr) = plot(variances([i j],dimCtr),repmat(distances(pairCtr,dimCtr),1,2),...
                    '.-','Color',plotColors(dimCtr,:));
                % label lines with strain letters
                if variances(i,dimCtr)<variances(j,dimCtr)
                    text(variances(i,dimCtr),distances(pairCtr,dimCtr),...
                        strains{i},'HorizontalAlignment','right')
                    text(variances(j,dimCtr),distances(pairCtr,dimCtr),...
                        strains{j},'HorizontalAlignment','left')
                else
                    text(variances(i,dimCtr),distances(pairCtr,dimCtr),...
                        strains{i},'HorizontalAlignment','left')
                    text(variances(j,dimCtr),distances(pairCtr,dimCtr),...
                        strains{j},'HorizontalAlignment','right')
                end
            end
        end
    end
    % annotate and save figure
    title(['N = ' num2str(N) ' worms'],'FontWeight','normal')
    box on
    xlabel('variance explained'), ylabel('difference in distributions (D_{JS}^{1/2})','Interpreter','tex')
    lg = legend(lineHandles,num2str((1:4)'));
    lg.Location = 'best';
    lg.Title.String = 'a_i';
    figName = ['figures/' num2str(N) 'worms_projections_Djs.eps'];
    jsFig.PaperUnits = 'centimeters';
    exportfig(jsFig,figName,exportOptions)
    system(['epstopdf ' figName]);
    system(['rm ' figName]);
    %             close(eigProjectionFig)
    %% plot triangles
    triFig = figure; 
    hold on
    for dimCtr = 1:nComponents
        [x, y] = triangle(distances(1,dimCtr),distances(2,dimCtr),distances(3,dimCtr));
        plot(x + 0.1*(dimCtr - 1),y,'.-','Color',plotColors(dimCtr,:),'LineWidth',2)
        % plot strain labels on vertices with size proportional to variance
        for strainCtr = 1:nStrains
           text(x(strainCtr) + 0.1*(dimCtr - 1),y(strainCtr),strains{strainCtr},...
               'HorizontalAlignment',textalignments{strainCtr},...
               'FontSize',20*sqrt(variances(strainCtr,dimCtr)),'FontWeight','bold')
        end
    end
    % annotate and save figure
    title(['N = ' num2str(N) ' worms'],'FontWeight','normal')
    ax = gca;
    ax.Visible = 'off';
    ax.DataAspectRatio = [1 1 1];
    lg = legend(num2str((1:4)'));
    lg.Location = 'EastOutside';
    lg.Title.String = 'a_i';
    figName = ['figures/' num2str(N) 'worms_projections_Djs_triangles.eps'];
    triFig.PaperUnits = 'centimeters';
    exportfig(triFig,figName,exportOptions)
    system(['epstopdf ' figName]);
    system(['rm ' figName]);
end
% tilefigs()