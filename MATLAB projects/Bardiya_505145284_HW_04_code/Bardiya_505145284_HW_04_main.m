%Ojas Bardiya
%UID: 505145284
%Homework_04

clc; clear all; close all;

%We first give an option to solve problem 1 or 2
%1 refers to the Split-and-Average and 2 refers to the Runge-Kutta Radioactivity problem
problem_chosen = input('Please enter 1 or 2 for the respective problem to be solved\n');
%call a switch statement to execute the problem chosen
switch(problem_chosen)
    case (1)
        %Set Initial conditions
        x0 = [0, 0, 1, 1];
        y0 = [0, 1, 0, 1];
        w = [1, 2, 1];
        
        %Set the conditions for the while-loop
        error = 1;
        count = 0;
        max_iterations = 25;
        
        %Plot the initial figure
        figure(1);
        plot(x0, y0, 'bo', 'Markersize', 10, 'MarkerFacecolor', 'b');
        xlim([0 max(x0)])
        ylim([0 max(y0)])
        axis equal
        hold on
        
        x = x0;
        y = y0;
        
        %do until maximum node displacement is less than given threshold
        while error > 10^-3 && count <= max_iterations
            %Split
            xs = splitPts(x);
            ys = splitPts(y);
            
            %Average
            xa = averagePts(xs, w);
            ya = averagePts(ys, w);
            
            %Calculate maximum error
            dx = xa - xs;
            dy = ya - ys;
            error = max(sqrt(dx.^2 + dy.^2));
            
            %update x and y
            x = xa;
            y = ya;
            
            %update the count
            count = count + 1;
        end
        
        %Plot the final values
        plot(x, y, 'ko', 'Markersize', 5, 'MarkerFacecolor', 'k');
        hold off
    case (2)
        %Set initial values 
        t_0 = 0;
        t_f = 15;
        %Carbon-15 half-life
        t_half = 2.45;
        y0 = 1;
        
        
        %Array for values of time-step
        dt = [1, 0.1, 0.01];
        
        fprintf('    dt       RK1       RK2       RK4\n');
        %Loop through the values of dt
        for i = 1:1:length(dt)
            %Set array for time values
            t = t_0:dt(i):t_f;
            %Calculate number of time steps
            nt = length(t);
            %Set arrays for different RK methods
            ye = zeros(1, nt);
            y1 = zeros(1, nt);
            y2 = zeros(1, nt);
            y4 = zeros(1, nt);
            
            %Set initial amount of carbon-15
            ye(1) = y0;
            y1(1) = y0;
            y2(1) = y0;
            y4(1) = y0;
            
            %Apply the different Runge-Kutta methods for each time-step
            for k = 1:1:nt-1
                ye(k+1) = ye(1)*(exp((-log(2)/t_half)*t(k+1)));
                %First Runge-Kutta
                y1(k+1) = advanceRK(y1(k), dt(i), 1);
                %Second Runge-Kutta
                y2(k+1) = advanceRK(y2(k), dt(i), 2);
                %Fourth Runge-Kutta
                y4(k+1) = advanceRK(y4(k), dt(i), 4);
            end
            
            %Calculate the mean error for each Runge-Kutta method
            fprintf(' %1.2f:  %.2e  %.2e  %.2e \n', dt(i), abs(mean(ye - y1)), abs(mean(ye - y2)), abs(mean(ye - y4)));
            %Plot the final values for different values of dt
            figure(i)
            plot(t, ye, 'k', 'LineWidth', 2);
            hold on
            plot(t, y1, 'r--', 'LineWidth', 3);
            plot(t, y2, 'b--', 'LineWidth', 3);
            plot(t, y4, 'g--', 'LineWidth', 3);
            xlabel('Time (s)');
            ylabel('Carbon-15');
            title(sprintf('Decay of Carbon 15 with a time-step of %1.2f',  dt(i)));
            legend({'Exact', 'RK-1', 'RK-2', 'RK-4'}, 'Location', 'northeast');
            
        end
        
    otherwise
        fprintf('Error: Please input 1 or 2 to choose which problem to solve.\n');
end