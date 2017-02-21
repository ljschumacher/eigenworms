function [ x,y ] = triangle(AB,AC,BC)
%TRIANGLE give coordinates of vertices given side lengths 
% the first coordinate is repeated as the fourth for easy plotting
x = zeros(4,1);
y = x;
x(2) = AB;
s = (AB + AC + BC)/2;
% calculate altitude in terms of sides, see https://en.wikipedia.org/wiki/Altitude_(triangle)#Altitude_in_terms_of_the_sides
h = 2*sqrt(s*(s - AB)*(s - AC)*(s - BC))/AB;
y(3) = h;
x(3) = AC*cos(asin(h/AC));
end

