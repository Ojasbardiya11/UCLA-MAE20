%Ojas Bardiya
%UID: 505145284
%Homework_05

clc; clear all; close all;

%We first give an option to solve problem 1 or 2
%1 refers to the Share-Birthday and 2 refers to the Random-Walk Radioactivity problem
problem_chosen = input('Please enter 1 or 2 for the respective problem to be solved\n');
%call a switch statement to execute the problem chosen
switch(problem_chosen)
    case (1)
        %Set number of trials
        n_trials = 1e4;
        %Initialize array for storing results of a trial
        N = zeros(n_trials, 1);
        
        %Monte-Carlo simulation
        for i = 1:1:n_trials
            %Set variable to determine if a match occurs
            is_match = 0;
            %Generate a random birthday
            group = ceil(rand*365);
            %Set variables for looping
            Itermax = 53;
            iter = 1;
            while is_match == 0 && iter <= Itermax
                iter = iter + 1;
                %Generate a new random birthday
                n_bday = ceil(rand*365);
                %Determine if a match exists
                for j = 1:1:length(group)
                    if abs(group(j) - n_bday) < 7 || abs(abs(group(j) - n_bday) - 365) < 7
                        is_match = 1;
                    end
                end
                %Add the new birthday to the group
                group = [group;n_bday];
            end
            %Enter the number of birthdays to find a match in the trial 
            N(i) = iter;
        end
        mid = median(N);
        fprintf('Median Number of People = %d\n', mid);
        %Create histogram
        figure(1)
        days = 2:1:53;
        histogram(N, days);
        grid on
        title('Distribution of trial results', 'Fontsize', 24);
    case (2)
        %Set Boundary conditions
        bounds = [5, -5, -5, 5];
        %Set number of trials
        n_trials = 5000;
        %Initialize array to hold number of moves per trial
        N_moves = zeros(n_trials, 1);
        %Set initial conditions
        xA0 = -5;
        yA0 = 0;
        xB0 = 5;
        yB0 = 0;
        %Start the Monte-Carlo simulation
        for i = 1:1:n_trials
            %Set initial conditions
            xAk = xA0;
            yAk = yA0;
            xBk = xB0;
            yBk = yB0;
            %set number of moves
            k = 0;
            %set variable to determine if collision has occurred
            has_collided = 0;
            while has_collided == 0 && k < 1000
                %Execute a random-walk
                [xAkp1, yAkp1] = RandWalk_2D(xAk, yAk, bounds);
                [xBkp1, yBkp1] = RandWalk_2D(xBk, yBk, bounds);
                %Commented out section is for the visualization
                %{
                %Create particle A on grid for step (k)
                xa_pos = [xAk - 0.5, xAk + 0.5, xAk + 0.5, xAk - 0.5];
                ya_pos = [yAk - 0.5, yAk - 0.5, yAk + 0.5, yAk + 0.5];
                %Create particle A on grid for step (k + 1)
                xakp1_pos = [xAkp1 - 0.5, xAkp1 + 0.5, xAkp1 + 0.5, xAkp1 - 0.5];
                yakp1_pos = [yAkp1 - 0.5, yAkp1 - 0.5, yAkp1 + 0.5, yAkp1 + 0.5];
                
                
                
                %Create particle A on grid for step (k)
                xb_pos = [xBk - 0.5, xBk + 0.5, xBk + 0.5, xBk - 0.5];
                yb_pos = [yBk - 0.5, yBk - 0.5, yBk + 0.5, yBk + 0.5];
                %Create particle A on grid for step (k + 1)
                xbkp1_pos = [xBkp1 - 0.5, xBkp1 + 0.5, xBkp1 + 0.5, xBkp1 - 0.5];
                ybkp1_pos = [yBkp1 - 0.5, yBkp1 - 0.5, yBkp1 + 0.5, yBkp1 + 0.5];
                
                %Create a grid for visualization for each trial
                %Particle A is red at step k and blue at step k+1
                %Particle B is yellow at step k and green at step k+1
                
                figure (i)
                hold on 
                set(gca,'xtick', -5.5:1:5.5)
                set(gca,'ytick', -5.5:1:5.5)
                grid on
                xlim([-5.5 5.5])
                ylim([-5.5 5.5])
                axis square
                fill(xa_pos, ya_pos, 'r')
                fill(xakp1_pos, yakp1_pos, 'b')
                fill(xb_pos, yb_pos, 'y')
                fill(xbkp1_pos, ybkp1_pos, 'g')
                title('2D Random-walk collision', 'FontSize', 24)
                set(gcf, 'Position', [30 350 850 450])
                set(gca, 'LineWidth', 3, 'FontSize', 20);
                hold off
                %}
                
                %Update the positions of both particles
                xAk = xAkp1;
                yAk = yAkp1;
                xBk = xBkp1;
                yBk = yBkp1;
                k = k + 1;
                %Determine if a collision has occurred
                if xAk == xBk && yAk == yBk
                    has_collided = 1;
                    N_moves(i) = k;
                end
                %If no collision, set the moves taken for a collision to the chosen limit
                if k == 1000
                    N_moves(i) = k;
                end
            end
        end
        mid = median(N_moves);
        fprintf('Median = %d\n', mid);
        
    otherwise
        fprintf('Error: Please input 1 or 2 to choose which problem to solve.\n');
end