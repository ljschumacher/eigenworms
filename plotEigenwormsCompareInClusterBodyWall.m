% combine eigenworm plots from saved eigenworms for different worm numbers
% plotting eigenworms in real space rather than in angles

close all
clear

% figure export options
exportOptions = struct('Color','rgb');

nEigenworms = 4;
strains = {'N2', 'npr1'};
nStrains = length(strains);
wormnums = {'HD','40'};
analysisTypes = {'loneWorms','inCluster','smallCluster'};

% for reference, eigenworms from Brown et al. 2013
reference = load('../Motif_Analysis/eigenWorms.mat');
% change the order of first and third masterworm to match ours
% reference.eigenWorms = reference.eigenWorms([3 2 1 4:end], :);

for strainCtr = 1:nStrains
    S = strains{strainCtr};
    for numCtr = 1:length(wormnums)
        N = wormnums{numCtr};
        eigWormFig = figure;
        for analysisType = analysisTypes
            % load eigenworm data
            file = rdir(['results/eigenData_' S '_' N '_bodywall_' analysisType{1} '.mat']);
            if ~isempty(file)
                % load eigenworm analysis result
                load(file.name,'eigenWorms')
                % plot first few eigenworms in real space
                [x, y] = angles2xy(eigenWorms(1:nEigenworms,:));
                for eigCtr = 1:nEigenworms
                    subplot(ceil(nEigenworms/2),2,eigCtr)
                    hold on
                    % check if we need to flip eigenworm upside-down
                    if eigenWorms(eigCtr,:)*reference.eigenWorms(eigCtr,:)'<0
                        eigenWorms(eigCtr,:) = -eigenWorms(eigCtr,:);
                    end
                    plot(x(eigCtr,:),y(eigCtr,:)-mean(y(eigCtr,:)),'LineWidth',2)
                end
            else
                display(['No data for strain=' S ', worms=' num2str(N)])
            end
        end
        % plot reference eigenworms from Brown et al. 2013
        [xref, yref] = angles2xy(reference.eigenWorms(1:nEigenworms,:));
        for eigCtr = 1:nEigenworms
            subplot(ceil(nEigenworms/2),2,eigCtr)
            plot(xref(eigCtr,:),yref(eigCtr,:)- mean(yref(eigCtr,:)),'k--')
            title(num2str(eigCtr),'FontWeight','normal')
            box on
            xlabel('x'), ylabel('y')
            xlim([0 48])
            ylim([-5 5])
        end
        % annotate and save figure
        legH = legend(eigWormFig.Children(1),analysisTypes,'Location','SouthEast');
        %     legH.Title.String = 'N worms';
        set(eigWormFig, 'name', ['Strain ' S ' ' N ' eigenworms'])
        figFileName = ['figures/eigenworms_compareInCluster_' S '_' N '.eps'];
        exportfig(eigWormFig,figFileName,exportOptions)
        system(['epstopdf ' figFileName]);
        system(['rm ' figFileName]);
    end
end
tilefigs()