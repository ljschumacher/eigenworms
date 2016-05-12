% plot mean eigen-projections in multiple dimensions to compare strains
close all
clear

files = dir('*.mat');

nStrains = size(files,1);

meanProj = cell(nStrains,6);
strains = cell(nStrains,1);

legendStrings = cell(nStrains,1);

for strainCtr = 1:nStrains
    load(files(strainCtr).name)
    legendStrings{strainCtr} = strrep(strrep(strrep(files(strainCtr).name,...
        '_',' '),'histogram',''),'.mat','');
    for eigCtr = 1:6
        meanProj{strainCtr,eigCtr} = ...
            worm.posture.eigenProjection(eigCtr).histogram.data.mean.abs;
    end
    strains{strainCtr} = strainCtr*ones(size(meanProj{strainCtr,1},1),1);
end

parallelcoords(cell2mat(meanProj),'Group',cell2mat(strains))
ylabel('\langle|eigen projection|\rangle')
xlabel('eigenworm number')
xlim([1 6])