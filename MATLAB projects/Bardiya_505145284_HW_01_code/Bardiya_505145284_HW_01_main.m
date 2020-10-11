%Ojas Bardiya
%UID: 505145284
%Homework_01

clc; clear all; close all;

%We first give an option to solve problem 1 or 2
%1 refers to the oblate spheroid and 2 refers to linear indexing
problem_chosen = input('Please enter 1 or 2 for the respective problem to be solved\n');

%call a switch statement to execute the problem chosen
switch(problem_chosen)
    case (1)
        %user-input
        r1 = input('Enter the equatorial radius:\n');
        r2 = input('Enter the polar radius:\n');
        %check if values entered are valid
        if (r2 >= r1)
            fprintf('Error:The equatorial radius must be larger than the polar radius!\n');
        elseif (r2 <= 0 || r1 <= 0)
            fprintf('Error:Both radii must have a non-negative real value!\n');
        else
            g = acos(r2/r1);
            %compute the surface area using the given formula
            S_area = 2*pi*(r1^2 + (r2^2/sin(g))*log(cos(g)/(1 - sin(g))));
            fprintf('%s%10e\n','The surface area is: ',S_area);
            %compute the approximation
            approx = 4*pi*((r1 + r2)/2)^2;
            fprintf('%s%10e\n','The approximation is: ',approx);
            %find the discrepancy
            discrepancy = S_area - approx;
            %output 10 digits in scienctific notation for the discrepancy
            fprintf('%s%10e\n','The discrepancy is: ',discrepancy);
        end
    case (2)
        %user-input
        is_error = 0;
        M = input('Enter the number of rows:\n');
        %make sure no. of rows is entered correctly
        if (M ~= floor(M))
            fprintf('Error: The number of rows must be an integer!\n');
            is_error = 1;
        elseif (M < 2)
            fprintf('Error: There must be more than 2 rows!\n');
            is_error = 1;
        end
        N = input('Enter the number of columns:\n');
        %make sure no. of columns is entered correctly
        if (N ~= floor(N))
            fprintf('Error: The number of columns must be an integer!\n');
            is_error = 1;
        elseif (N < 2)
            fprintf('Error: There must be more than 2 columns!\n');
            is_error = 1;
        end
        P = input('Enter the chosen cell:\n');
        %make sure the cell is entered correctly
        if (P ~= floor(P))
            fprintf('Error: The number of cells must be an integer!\n');
            is_error = 1;
        elseif (P < 1 || P > M*N)
            fprintf('Error: P must be in the range [1,M*N]\n');
            is_error = 1;
        end
        if (is_error == 1)
            fprintf('Error(s) have incurred previously, please check!\n')
        else
            %convert to integers
            P = int8(P);
            M = int8(M);
            N = int8(N);
            %list of all possible neighbours
            P_up = P - 1;
            P_down = P + 1;
            P_right = P + M;
            P_left =  P - M;
            %upper left diagnol
            P_lu = P - M - 1;
            %upper right diagnol
            P_ru = P + M - 1;
            %lower left diagnol
            P_ll = P - M + 1;
            %lower right diagnol
            P_rl = P + M + 1;
            
            %P is the top-left corner
            if (P == 1)
                neighbours = [P_down P_right P_rl];
            %P is the bottom-left corner
            elseif (P == M)
                neighbours = [P_up P_ru P_right];
            %P is the top-right corner    
            elseif (P == M*(N-1) + 1)
                neighbours = [P_left P_ll P_down];
            %P is the bottom-right corner
            elseif (P == M*N)
                neighbours = [P_lu P_left P_up];
            %P is on the left edge
            elseif (P > 1 && P < M)
                neighbours = [P_up P_down P_ru P_right P_rl];
            %P is on the right edge
            elseif (P > M*(N-1) + 1 && P < M*N)
                neighbours = [P_lu P_left P_ll P_up P_down];
            %P is on the bottom side
            elseif (mod(P,M) == 0 && P ~= M && P ~= N*M)
                neighbours = [P_lu P_left P_up P_ru P_right];
            %P is on the top side
            elseif (mod(P,M) == 1 && P ~= 1 && P ~= (M*(N-1) + 1))
                neighbours = [P_left P_ll P_down P_right P_rl];
            %P is an interior cell
            else
                neighbours = [P_lu P_left P_ll P_up P_down P_ru P_right P_rl];
            end
        end
        %display the results
        fprintf('%s%d\n', 'Cell ID:   ', P);
        fprintf('%s', 'Neighbors:');
        fprintf(' %d ', neighbours);
        fprintf('\n');
    otherwise
        fprintf('Error: Please input 1 or 2 to choose which problem to solve.\n')
end
    
        
            
            
            
        
