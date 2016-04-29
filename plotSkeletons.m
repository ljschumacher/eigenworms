close all
clear

% load some worm skeleta from tracking data and plot them
filename = '/data1/linus/Recordings/Results/Camille_151119/CSTCTest_Ch4_19112015_121940_skeletons.hdf5';

% read out the size of the dataset
[~,h5Dims] = H5S.get_simple_extent_dims(H5D.get_space(H5D.open(...
    H5F.open(filename),'/skeleton')));
% NOTE: the long dimension is the first one of h5dims, but the last one of
% the dataset

% select random frame numbers (so we don't get consecutive ones)
nPlots = 100;
frameIDs = ceil(rand(1,nPlots)*h5Dims(1));

% plot the skeleta
figure, hold on
% for spacing work plots
dx = 50;
dy = 50;
nx = 10;

for nn=1:length(frameIDs)
    ii = frameIDs(nn);
    % load one frame
    skeleton = h5read(filename,'/skeleton',[1,1,ii],[2,49,1]);
    % plot in on the next position on the 'grid'
    x0 = mean(skeleton(1,:));
    y0 = mean(skeleton(2,:));
    plot(skeleton(1,:) - x0 + mod(nn,nx)*dx, ...
        skeleton(2,:) - y0 + floor(nn/nx)*dy)
end