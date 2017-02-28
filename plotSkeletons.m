function [] = plotSkeletons(skelData,nFrames,color)

if nargin<3
    color = 'k';
end
% select random frame numbers (so we don't get consecutive ones)
nSkeleta = size(skelData,3);
skelData = skelData(:,:,randperm(nSkeleta,nFrames));

% for spacing work plots
dx = 200;
dy = 200;
nx = ceil(sqrt(nFrames));

hold on
for skelCtr=1:nFrames
    % plot in on the next position on the 'grid'
    x0 = mean(skelData(1,:,skelCtr));
    y0 = mean(skelData(2,:,skelCtr));
    plot(squeeze(skelData(1,:,skelCtr)) - x0 + mod(skelCtr-1,nx)*dx, ...
        squeeze(skelData(2,:,skelCtr)) - y0 + floor((skelCtr-1)/nx)*dy,...
        'Color',color)
end
axis equal