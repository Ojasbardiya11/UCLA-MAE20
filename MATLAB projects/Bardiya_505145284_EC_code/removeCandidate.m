function [votes] = removeCandidate(votes,losingCandidate)
%This function takes the votes for each candidate in the form of an array
%and outputs the same with the removal of the lowest-ranked candidate

dims = size(votes);
%find the number of voters
num_voters = dims(1);
%find the number of candidates
num_candidates = dims(2);
%create a new array for after the lowest-ranked candidate is eliminated
new_votes = zeros(num_voters, num_candidates - 1);

%Loop through the candidates
for i = 1:num_voters
    p = 0;
    %Loop through the voters
    for j = 1:num_candidates
        %If the vote is not for the losing candidate
        if votes(i, j) ~= losingCandidate
            %update the vote count for the candidate
            p = p + 1;
            new_votes(i, p) = votes(i,j);
        end
    end
end

%update the voters array
votes = new_votes;

end
        


