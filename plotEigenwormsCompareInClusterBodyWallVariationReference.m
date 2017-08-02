% combine eigenworm plots from saved eigenworms for different worm numbers
% plotting eigenworms in real space rather than in angles
% plotting variation of eigenworms

close all
clear

% figure export options
exportOptions = struct('Color','rgb');

nEigenworms = 4;
strains = {'N2', 'npr1'};
nStrains = length(strains);
wormnums = {'HD','40','1W'};
analysisTypes = {'loneWorms','inCluster','smallCluster','leaveCluster'};
nVariations = 5;
weights = linspace(-1,1,nVariations)';

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
            file = rdir(['results/eigenData_' S '_' N '_bodywall_' analysisType '.mat']);
            if ~isempty(file)
                % load eigenworm analysis result
                load(file.name,'eigenWorms','masterProjections')
                % plot first few eigenworms in real space
                for eigCtr = 1:nEigenworms
                    subplot(nEigenworms,length(analysisTypes),(eigCtr-1)*length(analysisTypes)+dataCtr), hold on
                    set(gca,'ColorOrder',parula(nVariations),'Color','none')
                    sigma2 = 2*std(masterProjections(:,eigCtr));
                    [x, y] = angles2xy(sigma2*weights*reference.eigenWorms(eigCtr,:));
                    % plot variation of reference eigenworms in real space, centered on y=0
                    plot((x-mean(x,2))',(y-mean(y,2))','LineWidth',2)
                    if eigCtr==1
                        title(analysisType,'FontWeight','normal')
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
        figFileName = ['figures/refeigenworms_variation_' S '_' N '.eps'];
        exportfig(eigWormFig,figFileName,exportOptions)
        system(['epstopdf ' figFileName]);
        system(['rm ' figFileName]);
    end
end
%tilefigs()
