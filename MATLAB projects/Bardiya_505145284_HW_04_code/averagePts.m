function [ xa ] = averagePts (xs, w)
%This function uses a 1*3 weight vector to determine the average of a
%series of points in the form of split array via a predetermined formula

%Check to see if the sum of weights is valid
if sum(w) == 0
    fprintf('Error: The sum of weights must be not be zero.\n');
else
    %Normalize the weight array
    w = w/sum(w);
    %Find the length of the split array
    n = length(xs);
    for k = 1:1:n
        %Determine the position of the left and right neighbours
        if k == 1
            head = n;
        else
            head = k - 1;
        end
        if k == n
            tail = 1;
        else
            tail = k + 1;
        end
        
        %Apply the formula for getting the weighted average
        xa(k) = w(1)*xs(head) + w(2)*xs(k) + w(3)*xs(tail);
    end
end


