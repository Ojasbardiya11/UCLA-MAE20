%Ojas Bardiya
%UID: 505145284
%Homework_02

clc; clear all; close all;

%We first give an option to solve problem 1 or 2
%1 refers to the 3-species and 2 refers to the pocket-change problem
problem_chosen = input('Please enter 1 or 2 for the respective problem to be solved\n');

%call a switch statement to execute the problem chosen
switch(problem_chosen)
    case (1)
        %Define the coefficients for the Lotka-Volterra Equation
        a = 0.75;
        b = 1.5;
        c = 0.5;
        d = 1.25;
        
        %Set the time-stepping parameters
        dt = 0.005; t_final = 12;
        
        %calculate the integer number of steps
        n_steps = ceil(t_final/dt);
        
        %initial parameters
        x_0 = 2; x_i = x_0;
        y_0 = 2.49; y_i = y_0;
        z_0 = 1.5; z_i = z_0;
        
        %ouput the initial conditions
        fprintf('Time    X    Y    Z\n');
        fprintf(' %1.1f %1.2f %1.2f %1.2f\n',i*dt,x_i,y_i,z_i);
        
        %Begin the for-loop
        for i = 1:1:n_steps
            
            %Forward Euler updating method
            x_ip1 = x_i + dt*(a*x_i*(1 - x_i/20) - b*x_i*y_i - c*x_i*z_i);
            y_ip1 = y_i + dt*(y_i*(1 - y_i/25) - a*x_i*y_i - d*y_i*z_i);
            z_ip1 = z_i + dt*(b*z_i*(1 - z_i/30) - x_i*z_i - y_i*z_i);
            
            %update values
            x_i = x_ip1;
            y_i = y_ip1;
            z_i = z_ip1;
            
            %check if values should be printed
            if mod(i*dt,0.5) == 0
                fprintf(' %1.1f %1.2f %1.2f %1.2f\n',i*dt,abs(x_i),abs(y_i),abs(z_i));
            end
        end
    case (2)
        %Set the total count of coins
        coin_count = 0;
        for Moneycount = 0:99
            %Set total count for each denomination
            %quarters
            Q = 0;
            %dimes
            D = 0;
            %nickels
            N = 0;
            %pennies
            P = 0;
            %Set remaining amount of money 
            Remaining_amount = Moneycount;
            %use a greedy approach to optimize the coin selection
            while Remaining_amount > 0
                if Remaining_amount >= 25
                    Q = Q + 1;
                    Remaining_amount = Remaining_amount - 25;
                elseif Remaining_amount >= 10
                    D = D + 1;
                    Remaining_amount = Remaining_amount - 10;
                elseif Remaining_amount >= 5
                    N = N + 1;
                    Remaining_amount = Remaining_amount - 5;
                elseif Remaining_amount >= 1
                    P = P + 1;
                    Remaining_amount = Remaining_amount - 1;
                end
            end
            %update the coin count
            coin_count = coin_count + Q + D + N + P;
        end
        %calculate the average 
        Average_coins = coin_count/100;
        fprintf('Average Number of Coins = %.2f\n', Average_coins);
    otherwise
        fprintf('Error: Please input 1 or 2 to choose which problem to solve.\n');
end