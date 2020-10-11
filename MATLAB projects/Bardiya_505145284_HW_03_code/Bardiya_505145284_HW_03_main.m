%Ojas Bardiya
%UID: 505145284
%Homework_03

clc; clear all; close all;

%We first give an option to solve problem 1 or 2
%1 refers to the pendulum and 2 refers to the DNA-analysis problem
problem_chosen = input('Please enter 1 or 2 for the respective problem to be solved\n');
%call a switch statement to execute the problem chosen
switch(problem_chosen)
    case (1)
        %Setting up initial values
        theta_0 = pi/3;
        g = 9.81;
        L = 1;
         
        %Setting up time-step values
        t_0 = 0;
        t_f = 20;
        dt = 0.005;
        nt = (t_f - t_0)/dt;
        t = linspace(t_0, t_f, nt);
         
        %Initialize the arrays
        theta = zeros(1, nt); theta(1) = theta_0;
        omega = zeros(1, nt);
        alpha = zeros(1,nt);
        
        %Choose the desired method 
        %1 for Forward Euler and 2 for semi-implicit Euler
        method = input('Enter the approximation method:\n');
        switch(method)
            case (1)
                %Forward Euler Iteration
                for k = 1:1:nt-1
                    omega(k+1) = -g/L*sin(theta(k))*dt + omega(k);
                    theta(k+1) = omega(k)*dt + theta(k);
                    alpha(k+1) = (omega(k+1) - omega(k))/dt;
                end
            case (2)
                %Semi-implicit Euler Iteration
                for k = 1:1:nt-1
                    omega(k+1) = -g/L*sin(theta(k))*dt + omega(k);
                    theta(k+1) = omega(k+1)*dt + theta(k);
                    alpha(k+1) = (omega(k+1) - omega(k))/dt;
                end
            otherwise
                fprintf('Error: Please select 1 or 2 to choose a particular method!');
        end
        %Calculate the energy per unit mass
        E = g*L*(1 - cos(theta)) + 0.5*(L.*omega).^2;
        
        %Plots
        %Plot 1 - Kinematic vectors vs time
        figure(1)
        plot(t, theta, t, omega, t, alpha);
        xlabel('Time (s)');
        ylabel('Kinematic vectors');
        title('Plotting Kinematic vectors over time for a pendulum');
        legend('Angular Position', 'Angular Speed', 'Accelaration');
        
        %Plot 2 - Energy per unit mass vs time
        figure(2)
        plot(t, E);
        xlabel('Time (s)');
        ylabel('Energy per unit mass (J/kg)');
        title('Plotting energy per unit mass over time for a pendulum');
        
    case (2)
        %DNA-analysis
        file = load('chr1_sect.mat');
        dna = file.dna;
        
        %Set count for each stop codon
        taa = 0;
        tag = 0;
        tga = 0;
        
        %Initialize the startpoint as zero
        start = 0;
        
        %get length of the chromosome
        numbases = length(dna);
        
        %Create an array to hold the value of the segment lengths
        length_values = zeros(1, ceil(numbases/6));
        cur_index = 0;
        
        %Loop through the bases in the chromosome
        for k = 1:3:numbases - 2
            if start == 0
                if dna(k) == 1 && dna(k+1) == 4 && dna(k+2) == 3
                    start = k;
                end
            
            
            else
                %check for TGA
                if dna(k) == 4 && dna(k+1) == 3 && dna(k+2) == 1 
                    tga = tga + 1;
                    seg_length = k - start + 3;
                    cur_index = cur_index + 1;
                    length_values(cur_index) = seg_length;
                    start = 0;
                %check for TAA
                elseif dna(k) == 4 && dna(k+1) == 1 && dna(k+2) == 1 
                    taa = taa + 1;
                    seg_length = k - start + 3;
                    cur_index = cur_index + 1;
                    length_values(cur_index) = seg_length;
                    start = 0;
                %check for TAG 
                elseif dna(k) == 4 && dna(k+1) == 1 && dna(k+2) == 3 
                    tag = tag + 1;
                    seg_length = k - start + 3;
                    cur_index = cur_index + 1;
                    length_values(cur_index) = seg_length;
                    start = 0;
                end
            end
        end
        
        %Calculate the total segments
        ts = cur_index;
        %Calculate the average length
        avg_length = mean(length_values(1:ts));
        %Calculate the minimum length
        min_length = min(length_values(1:ts));
        %Calculate the maximum length
        max_length = max(length_values(1:ts));
        
        %Print output accordingly
        fprintf('Total protein coding segments: %d\n', ts);
        fprintf('Average length: %.2f\n', avg_length);
        fprintf('Maximum length: %d\n', max_length);
        fprintf('Minimum length: %d\n', min_length);
        
        
    otherwise
        fprintf('Error: Please input 1 or 2 to choose which problem to solve.\n');
end
        
                
            
            
                
              
        
        
