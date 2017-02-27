function [ x, y ] = angles2xy( angles )
% converts vector of (eigenworm) tangent angles back into xy coordinates
% assumes the longer dimension is the one to sum
if size(angles,2) >= size(angles,1)
    dimsum = 2;
    numSamples = size(angles,1);
else
    dimsum = 1;
    numSamples = size(angles,2);
end
x = [zeros(numSamples,1) cumsum(cos(angles),dimsum)];
y = [zeros(numSamples,1) cumsum(sin(angles),dimsum)];
end

