% combine eigenworm plots from saved eigenworms for different worm numbers
% plotting eigenworms in real space rather than in angles
% plotting variation of eigenworms

close all
clear

% specify how to phase-restrict
phase = 'joining'; % 'fullMovie', 'joining', or 'sweeping'.

% figure export options
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',1,...
    'Renderer','painters');

nEigenworms = 4;
strains = {'N2', 'npr1'};
nStrains = length(strains);
wormnums = {'40'}%{'HD','40','1W'};
if ~strcmp(phase, 'fullMovie')
    wormnums = {'40'};
end
analysisTypes = {'loneWorms','inCluster'}%,'smallCluster','leaveCluster'};
analysisLabels = {'lone worms', 'in cluster'};
quants = [0.2, 0.4, 0.6, 0.8];
nVariations = length(quants);

% for reference, eigenworms from Brown et al. 2013
reference = load('eigenWorms.mat');

% % change the order of first and third masterworm
% reference.eigenWorms = reference.eigenWorms([3 2 1 4:end], :);

%% also use this to plot variation on projections of master eigenworms
for strainCtr = 1:nStrains
    S = strains{strainCtr};
    for numCtr = 1:length(wormnums)
        N = wormnums{numCtr};
        eigWormFig = figure;
        for dataCtr = 1:length(analysisTypes)
            analysisType = analysisTypes{dataCtr};
            % load eigenworm data
            file = rdir(['results/eigenData_' S '_' N '_bodywall_' analysisType '_' phase '.mat']);
            if ~isempty(file)
                % load eigenworm analysis result
                load(file.name,'eigenWorms','masterProjections')
                % plot first few eigenworms in real space
                for eigCtr = 1:nEigenworms
                    subplot(nEigenworms,length(analysisTypes),(eigCtr-1)*length(analysisTypes)+dataCtr), hold on
                    set(gca,'ColorOrder',parula(nVariations),'Color','none')
                    weights = quantile(abs(masterProjections(:,eigCtr)),quants);
                    [x, y] = angles2xy(weights'*reference.eigenWorms(eigCtr,:));
                    % plot variation of reference eigenworms in real space, centered on y=0
                    plot((x-mean(x,2))',(y-mean(y,2))','LineWidth',2)
                    if eigCtr==1
                        title(analysisLabels{dataCtr},'FontWeight','normal')
                    end
                    if numCtr==1
                        ylabel(num2str(eigCtr))
                    end
                    xlim([-24 24])
                    axis equal
                    set(gca,'XTick',[],'YTick',[])
                end
            else
                display(['No data for strain=' S ', worms=' num2str(N) ', ' analysisType])
            end
        end
        % annotate and save figure
        set(eigWormFig, 'name', ['Strain ' S ' ' N ' reference eigenworms'])
        figFileName = ['figures/refeigenworms_variation_' S '_' N '_' phase '.eps'];
        exportfig(eigWormFig,figFileName,exportOptions)
        system(['epstopdf ' figFileName]);
        system(['rm ' figFileName]);
    end
end
%tilefigs()
