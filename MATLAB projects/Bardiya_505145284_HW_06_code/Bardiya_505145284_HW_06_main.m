%Ojas Bardiya
%UID: 505145284
%Homework_06

clc; clear all; close all;

%We first give an option to solve problem 1 or 2
%1 refers to the Game of Life and 2 refers to the Euler-Bernoulli Beam Bending problem
problem_chosen = input('Please enter 1 or 2 for the respective problem to be solved\n');
%call a switch statement to execute the problem chosen
switch(problem_chosen)
    case (1)
        %Set the grid size
        %set no. of columns
        x_grid = 200;
        %set no. of rows
        y_grid = 150;
        
        %Initialize LifeGrid
        LifeGrid = rand(y_grid, x_grid);
        for i = 1:1:y_grid
            for j = 1:1:x_grid
                if LifeGrid(i, j) <= 0.1
                    LifeGrid(i, j) = 1;
                else
                    LifeGrid(i, j) = 0;
                end
            end
        end
        
        %Visualize the initial condition
        figure(1)
        imagesc(LifeGrid)
        title('Game of Life at n = 0', 'FontSize', 24)
        set(gcf, 'Position', [30 350 850 450])
        set(gca, 'LineWidth', 3, 'FontSize', 20);
        axis equal
        drawnow
        
        %Set the number of generations
        n_gens = 300;
        %Intialize array to hold count of living cells in each generation
        count = zeros(1, n_gens);
        
        
        %Iterate the Game of Life via time-steps
        iter = 0;
        while iter < n_gens
            pause(0.01)
            iter = iter + 1;
            
            %Create a new grid to hold values for the next time-step
            NewGrid = LifeGrid;
            %Loop through the rows
            for j = 1:1:y_grid
                %check for boundary conditions
                if j == 1
                    jm1 = y_grid;
                else
                    jm1 = j - 1;
                end
                if j == y_grid
                    jp1 = 1;
                else
                    jp1 = j + 1;
                end
                %Loop through the columns
                for i = 1:1:x_grid
                    %check boundary conditions
                    if i == 1
                        im1 = x_grid;
                    else
                        im1 = i - 1;
                    end
                    if i == x_grid
                        ip1 = 1;
                    else
                        ip1 = i + 1;
                    end 
                    %Get the sum of neighbours for the current cell
                    l_neighbours = LifeGrid(jm1, im1) + LifeGrid(jm1, i) + LifeGrid(jm1,ip1)+...
                                   LifeGrid(j, im1) + LifeGrid(j, ip1) + ...
                                   LifeGrid(jp1, im1) + LifeGrid(jp1, i) + LifeGrid(jp1, ip1);
                    %if the cell is alive
                    if LifeGrid(j,i) == 1
                        %Rule 1: A living cell survives if it has 2 or 3
                        %living neighbours
                        if l_neighbours == 2 || l_neighbours == 3
                            NewGrid(j,i) = 1;
                        %Rule 2: A living cell having < 2 or > 3 neighbours
                        %does not survive into the next generation
                        else
                            NewGrid(j,i) = 0;
                        end
                    %the cell is dead
                    else
                        %Rule 3: A dead cell becomes alive in the next generation if it has
                        %exactly 3 neighbours
                        if l_neighbours == 3
                            NewGrid(j,i) = 1;
                        end
                    end
                end
            end
            
            %Plot the current Life Grid
            imagesc(NewGrid)
            title(['The Game of Life at n = ' num2str(iter) '.'],'FontSize', 24);
            set(gcf, 'Position', [30 350 850 450])
            set(gca, 'LineWidth', 3, 'FontSize', 20);
            axis equal
            drawnow
            
            %Update the current LifeGrid Matrix
            LifeGrid = NewGrid;
            
            %sum the total number of living cells
            for m = 1:1:y_grid
                for n = 1:1:x_grid
                    count(iter) = count(iter) + LifeGrid(m, n);
                end
            end
        end
        
        %Plot no. of living cells in each generation
        figure(2)
        plot(linspace(0,n_gens,n_gens),count, '-');
        xlabel('Current Generation');
        ylabel('No. of Living Cells');
        
        
    case (2)
        %Set initial conditions
        
        %outer radius
        R = 0.013;
        %inner radius
        r = 0.011;
        %beam length
        L = 1;
        %Young's modulus
        E = 70e9;
        %distance
        d = 0.75;
        %Applied force 
        P = 2000;
        %Moment of inertia 
        I = (pi/4)*(R^4 - r^4);
        %no. of nodes
        N = 20;
        %Initialize the node coordinates
        x = zeros(N, 1);
        dx = L/(N - 1);
        for i = 2:1:N
            x(i) = x(i - 1) + dx;
        end
        %Set up a matrix to discretize the second-order differential
        %equation
        A = zeros(N,N);
        A(1, 1) = 1;
        A(N, N) = 1;
        for i = 2:1:N-1
            A(i, i - 1) = 1;
            A(i, i) = -2;
            A(i, i + 1) = 1;
        end
        
        %Calculate the bending moment 
        M = zeros(N, 1);
        for i = 1:N
            if x(i) <= d
                M(i) = -P*(L - d)*x(i)/L;
            else
                M(i) = -P*d*(L - x(i))/L;
            end
        end
        
        %solution vector
        b = dx^2*M/(E*I);
        %solve for y in the matrix system Ay = b
        y = A\b;
        %plot displacement vs x-position
        plot(x, y, '-o');
        xlabel('x-position (m)');
        ylabel('Displacement (m)');
        
        %Calculation for error
        c = min(d,L - d);
        %Find maximum displacement according to the calculated result
        y_max1 = max(y);
        %Find maximum displacement according to the theoretical solution
        y_max2 = (P*c*(L^2 - c^2)^1.5)/(9*sqrt(3)*E*I*L);
        %Determine error
        error = (abs(y_max2 - y_max1))/y_max2;
        
        fprintf('The error observed is %.2f percent\n', error*100);
        
        
        
    otherwise
        fprintf('Error: Please input 1 or 2 to choose which problem to solve.\n');
end