% calculate JS metric between distributions of eigencomponents to
% distinguish different strains

% issues/ to-do
% - is not yet fully flexible for different numbers of components, eg
% legendstring
close all
clear

% figure export options
exportOptions = struct('Color','rgb','Width',15,'FontMode','scaled','FontSizeMin',1);

nComponents = 4;
nComponentPairs = nComponents*(nComponents - 1)/2;
binLimits = [-5 5];
binWidth = 0.1;
nBins = length(min(binLimits):binWidth:max(binLimits));
strains = {'N2', 'HW', 'NP'};
textalignments = {'right','left','center'};
lineHandles = NaN(nComponentPairs,1);
legendstring = ['1-2';'1-3';'1-4';'2-3';'2-4';'3-4'];
nStrains = length(strains);
plotColors = lines(nComponentPairs);
for N = [40]
    distances = NaN(nStrains*(nStrains - 1)/2,nComponentPairs);
    distributions =  NaN(nStrains,nComponentPairs,nBins - 1,nBins - 1);
    variances = NaN(nStrains,nComponentPairs);
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
            masterVarExplained = masterVarExplained([3 2 1 4:nComponents]);
            % normalise to unit variance
            masterProjections = zscore(masterProjections);
            % loop through projected components
            for ii = 1:nComponents
                for jj = (ii + 1):nComponents
                % calculate 2D histograms
                dim2dCtr = (ii - 1)*(nComponents - ii/2) + jj - ii;
                distributions(strainCtr,dim2dCtr,:,:) = histcounts2(masterProjections(:,ii),masterProjections(:,jj),...
                    'BinWidth',binWidth,'XBinLimits',binLimits,'YBinLimits',binLimits,'Normalization','Probability');
                variances(strainCtr,dim2dCtr) = masterVarExplained(ii) + masterVarExplained(jj);
                end
            end
            clear masterProjections
        else
            display(['No data for strain=' S ', worms=' num2str(N)])
        end
    end
    % for each component calculate js-metric between all pairs of strains
    for dim2dCtr = 1:nComponentPairs
        for i = 1:nStrains
            for j = (i + 1):nStrains
                pairCtr = (i - 1)*(nStrains - i/2) + j - i;
                distances(pairCtr,dim2dCtr) = ...
                    sqrt(jsdiv(distributions(i,dim2dCtr,:,:),distributions(j,dim2dCtr,:,:)));
                % plot line at height=distance and length = difference in
                % variance explained
                lineHandles(dim2dCtr) = plot(variances([i j],dim2dCtr),repmat(distances(pairCtr,dim2dCtr),1,2),...
                    '.-','Color',plotColors(dim2dCtr,:));
                % label lines with strain letters
                if variances(i,dim2dCtr)<variances(j,dim2dCtr)
                    text(variances(i,dim2dCtr),distances(pairCtr,dim2dCtr),...
                        strains{i},'HorizontalAlignment','right')
                    text(variances(j,dim2dCtr),distances(pairCtr,dim2dCtr),...
                        strains{j},'HorizontalAlignment','left')
                else
                    text(variances(i,dim2dCtr),distances(pairCtr,dim2dCtr),...
                        strains{i},'HorizontalAlignment','left')
                    text(variances(j,dim2dCtr),distances(pairCtr,dim2dCtr),...
                        strains{j},'HorizontalAlignment','right')
                end
            end
        end
    end
    % annotate and save figure
    title(['N = ' num2str(N) ' worms'],'FontWeight','normal')
    box on
    xlabel('variance explained'), ylabel('difference in 2d distributions (D_{JS}^{1/2})','Interpreter','tex')
    lg = legend(lineHandles,legendstring);
    lg.Location = 'best';
    lg.Title.String = 'a_i-a_j';
    figName = ['figures/' num2str(N) 'worms_projections_Djs_2D.eps'];
    jsFig.PaperUnits = 'centimeters';
    exportfig(jsFig,figName,exportOptions)
    system(['epstopdf ' figName]);
    system(['rm ' figName]);
    %             close(eigProjectionFig)
    %% plot triangles
    triFig = figure; 
    hold on
    for dim2dCtr = 1:nComponentPairs
        [x, y] = triangle(distances(1,dim2dCtr),distances(2,dim2dCtr),distances(3,dim2dCtr));
        plot(x + 0.6*(dim2dCtr - 1),y,'.-','Color',plotColors(dim2dCtr,:),'LineWidth',2)
        % plot strain labels on vertices with size proportional to variance
        for strainCtr = 1:nStrains
           text(x(strainCtr) + 0.6*(dim2dCtr - 1),y(strainCtr),strains{strainCtr},...
               'HorizontalAlignment',textalignments{strainCtr},...
               'FontSize',20*sqrt(variances(strainCtr,dim2dCtr)),'FontWeight','bold')
        end
    end
    % annotate and save figure
    title(['N = ' num2str(N) ' worms'],'FontWeight','normal')
    ax = gca;
    ax.Visible = 'off';
    ax.DataAspectRatio = [1 1 1];
    lg = legend(legendstring);
    lg.Location = 'EastOutside';
    lg.Title.String = 'a_i-a_j';
    figName = ['figures/' num2str(N) 'worms_projections_Djs_2D_triangles.eps'];
    triFig.PaperUnits = 'centimeters';
    exportfig(triFig,figName,exportOptions)
    system(['epstopdf ' figName]);
    system(['rm ' figName]);
end
% tilefigs()