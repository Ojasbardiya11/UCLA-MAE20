%Ranked-choice voting
%Ojas Bardiya
%505145284

clc; clear all; close all;

%load the file data
file = load('votes1.mat');
votes = file.votes;

%Set up variables
winner_flag = 0;
cur_round = 0;

dims = size(votes);
%find the number of voters
num_voters = dims(1);
%find the number of candidates
num_candidates = dims(2);

fprintf('\t\t\t');
for k = 1:num_candidates
    fprintf('%d        ', k);
end
fprintf('\n');


while winner_flag ~= 1

    %Tally the votes for each candidate
    cur_round = cur_round + 1;
    total_count = zeros(1, num_candidates);
    for i = 1:num_voters
        k = votes(i, 1);
        total_count(k) = total_count(k) + 1;
    end
    
    %Determine if a majority vote has been established
    M = max(total_count)/num_voters;
    if M > 0.5
        winner_flag = 1;
    end
    
    %Determine the winning and losing candidate
    [total_sorted, candidate_num] = sort(total_count);
    winning_can = candidate_num(num_candidates);
    losing_can = candidate_num(cur_round);
    
    %remove the least ranked candidate
    votes=removeCandidate(votes, losing_can);

    %Print the current round results
    fprintf('Round %d Totals:        ', cur_round);
    for p = 1:num_candidates
        fprintf('%4d     ', total_count(p))
    end
    fprintf('\n');
    
end

%Display the winner
fprintf('Winning Candidate:  %d\n', winning_can);


%Extra credit part (c)
%{
If we break the tie in round 4 in favour of Candidate 2, we see the final
tally for Candidate 2 is 2353 and Candidate 3 is 2647. 
If we break the tie in round 4 in favour of Candidate 4, we see that the final tally for
Candidate 4 is 2400 and Candidate 3 is 2600. 
Thus, in both cases we see that 3 is the winning Candidate but wins by a larger margin if we choose candidate 2
as the winner of tie break in place of Candidate 4. 
Thus, while it does not have an effect on the final outcome, it does affect the final tally of
votes; in other circumstances, there may be a case where choosing a particular candidate 
in a tie break results in a different winner.
%}



