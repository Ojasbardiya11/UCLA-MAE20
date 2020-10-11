function [ xs ] = splitPts ( x )
%This function uses the midpoint formula to calculate the midpoint of every
%two consecutive elements in an array and creates a new array which
%contains both the original points and newly inserted midpoints under the
%assumption that the first and last element 'wrap around' each other, i.e.,
%the array is a circle
n = length(x);
for i = 1:1:n
    %Populate new points with the previous ones for every other index
    xs(2*i - 1) = x(i);
    %If the index is the last,assume it cycles to the first again
    if i == n
        xs(2*i) = (x(i) + x(1))/2;
    %Populate remaining cells at the midpoint of every pair of adjacent cells    
    else
        xs(2*i) = (x(i) + x(i + 1))/2;
    end
end



