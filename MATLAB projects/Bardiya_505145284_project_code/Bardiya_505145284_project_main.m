%Ojas Bardiya
%UID: 505145284
%Final Project - Mass-Spring Damper set up

clc; clear all; close all;

%Setting up initial values
%mass
m = [3 4 5 8 10 6 20 12];
%Spring constant
k = [200 50 125 25 300 100 80 75];
%Damping ratio
c = [2 45 50 35 10 4 80 65];

%Time parameters
t0 = 0;
tf = 10;
dt = 1/300;
t = t0:dt:tf;
nt = length(t);

%number of cases
%Outputs a 8*nt matrix containing values for each test case
%values 6-8 for m,k and c are user-defined
%6 corresponds to under-damped, 7 corresponds to critically-damped, 8
%corresponds to over-damped vibration.
n_exp = 8;

%Natural Frequency
w_n = sqrt(k./m);
%Damping ratio
xi = c./(2.*m.*w_n);

%Initial values
x0 = 1;
a0 = 2.5;

%Free Vibration - 1st Order Runge-Kutta array Initialization
x_RK1_free = zeros(n_exp, nt);
v_RK1_free = zeros(n_exp, nt);

%Free Vibration - 2nd Order Runge Kutta array Initialization
x_RK2_free = zeros(n_exp, nt);
v_RK2_free = zeros(n_exp, nt);

%Free Vibration - 4th Order Runge Kutta array Initialization
x_RK4_free = zeros(n_exp, nt);
v_RK4_free = zeros(n_exp, nt);


for i = 1:1:n_exp
    
    %Initial condition
    x_RK1_free(i, 1) = x0;
    x_RK2_free(i, 1) = x0;
    x_RK4_free(i, 1) = x0;
    
    for j = 1:1:nt-1
        %Free vibration - force is 0
        f = [0 0 0];
        
        %Find the [position, velocity] for the next iteration
        
        %Runge-Kutta 1
        temp = VibrationPosition([x_RK1_free(i, j) v_RK1_free(i, j)], m(i), k(i), c(i), f, dt, 1);
        %Runge-Kutta 2
        temp2 = VibrationPosition([x_RK2_free(i, j) v_RK2_free(i, j)], m(i), k(i), c(i), f, dt, 2);
        %Runge-Kutta 4
        temp4 = VibrationPosition([x_RK4_free(i, j) v_RK4_free(i, j)], m(i), k(i), c(i), f, dt, 4);
        
        %Update the values
        x_RK1_free(i, j+1) = temp(1);
        v_RK1_free(i, j+1) = temp(2);
        
        x_RK2_free(i, j+1) = temp2(1);
        v_RK2_free(i, j+1) = temp2(2);
        
        x_RK4_free(i, j+1) = temp4(1);
        v_RK4_free(i, j+1) = temp4(2);
    end
end


%Forced Vibration - 1st Order Runge-Kutta array Initialization
x_RK1_force = zeros(n_exp, nt);
v_RK1_force = zeros(n_exp, nt);

%Forced Vibration - 2nd Order Runge-Kutta array Initialization
x_RK2_force = zeros(n_exp, nt);
v_RK2_force = zeros(n_exp, nt);

%Forced Vibration - 4th Order Runge-Kutta array Initialization
x_RK4_force = zeros(n_exp, nt);
v_RK4_force = zeros(n_exp, nt);

for i = 1:1:n_exp
    
    %Initial condition
    x_RK1_force(i, 1) = x0;
    x_RK2_force(i, 1) = x0;
    x_RK4_force(i, 1) = x0;
    
    for j = 1:1:nt-1
        %Force vector
        f = [a0*sin(t(j)/(2*pi))/m(i) a0*sin((t(j) + 0.5*dt)/(2*pi))/m(i) a0*sin((t(j) + dt)/(2*pi))/m(i)];
        
        %Find the [position, velocity] for the next iteration
        %Runge-Kutta 1
        temp = VibrationPosition([x_RK1_force(i, j) v_RK1_force(i, j)], m(i), k(i), c(i), f, dt, 1);
        %Runge-Kutta 2
        temp2 = VibrationPosition([x_RK2_force(i, j) v_RK2_force(i, j)], m(i), k(i), c(i), f, dt, 2);
        %Runge-Kutta 4
        temp4 = VibrationPosition([x_RK4_force(i, j) v_RK4_force(i, j)], m(i), k(i), c(i), f, dt, 4);
        
        %Update the values
        x_RK1_force(i, j+1) = temp(1);
        v_RK1_force(i, j+1) = temp(2);
        
        x_RK2_force(i, j+1) = temp2(1);
        v_RK2_force(i, j+1) = temp2(2);
        
        x_RK4_force(i, j+1) = temp4(1);
        v_RK4_force(i, j+1) = temp4(2);
    end
end

%Plotting the graphs
for i = 1:1:n_exp
    figure(i)
    hold on
    set(gcf, 'Position', [15 50 1350 775])
    subplot(1,2,1)
    hold on
    grid on
        plot(t, x_RK1_free(i,:), 'ro')
        plot(t, x_RK2_free(i,:), 'g--', 'LineWidth', 3)
        plot(t, x_RK4_free(i,:), 'b--', 'LineWidth', 4)
        set(gca, 'LineWidth', 3, 'FontSize', 20)
        title('Homogenous Response')
        legend({'Forward Euler', 'RK-2', 'RK-4'}, 'Location', 'northeast')
    subplot(1,2,2)
    hold on
    grid on
        plot(t, x_RK1_force(i,:), 'ro')
        plot(t, x_RK2_force(i,:), 'g--', 'LineWidth', 3)
        plot(t, x_RK4_force(i,:), 'b--', 'LineWidth', 4)
        set(gca, 'LineWidth', 3, 'FontSize', 20)
        legend({'Forward Euler', 'RK-2', 'RK-4'}, 'Location', 'northeast')
        title('Inhomogenous Response')
end


%Frequency response
%Determine normalized frequency
lambda = -logspace(-1,1,500);
n_lambda = length(lambda);
%Intialize arrays absolute Gain and phase shift
Gain = zeros(n_exp, n_lambda);
PhiG = zeros(n_exp, n_lambda);

figure(9)
set(gcf, 'Position', [15 50 1350 775])
for run = 1:1:n_exp
    
    for j = 1:1:n_lambda
        %Determine values for absolute Gain and phase shift
        Gain(run, j) = 1/sqrt((1-lambda(j)^2)^2 + (2*xi(run)*lambda(j))^2);
        PhiG(run, j) = -atan(2*xi(run)*lambda(j)/(1-lambda(j)^2));
    end
    
    %Plot the relevant graphs
    subplot(2,1,1)
        semilogx(lambda,Gain(run,:), 'Linewidth',3)
        hold on
        grid on
        set(gca, 'LineWidth', 3, 'FontSize', 20)
    subplot(2,1,2)
        semilogx(lambda,PhiG(run,:), 'Linewidth',3);
        hold on
        grid on
        set(gca, 'LineWidth', 3, 'FontSize', 20)
    
   %Create a dictionary with the different values of run as keys 
   legendInfo{run} = ['Run ' num2str(run)];
end
legend(legendInfo)


%Animation - Underdamped case
%define filename for video
video_filename = 'Underdamped'

%create video file
vidfile = VideoWriter(video_filename, 'MPEG-4');
vidfile.FrameRate = 30;
open(vidfile);

%Initialize frame
nv = ceil(length(t)/10);
t_out = zeros(nv, 1);
x_hom = zeros(nv, 1);
x_inh = zeros(nv, 1);
count = 0;
sl = 0.25;


figure(10)
set(gcf,'Position',[15 50 1350 775])

for i = 1:1:length(t)
    %Plot a frame every 10s
    if mod(i-1,10) == 0
        count = count + 1;
        t_out(count) = t(i);
        %Trial 6 corresponds to user-defined underdamped case
        x_hom(count) = x_RK1_free(6,i);
        x_inh(count) = x_RK1_force(6,i);
        
        %Define the axes in the clockwise direction
        y1 = [x_hom(count)+sl x_hom(count)-sl x_hom(count)-sl x_hom(count)+sl];
        y2 = [x_inh(count)+sl x_inh(count)-sl x_inh(count)-sl x_inh(count)+sl];
        x = [-sl -sl sl sl];
        
        subplot(1,2,1)
            fill(x, y1, 'b')
            xlim([-1-sl 1+sl])
            ylim([-1-sl 1+sl])
            xlabel('X location [m]')
            ylabel('Y location [m]')
            title('Homogenous response')
            axis square
        subplot(1,2,2)
            fill(x, y2, 'b')
            xlim([-1-sl 1+sl])
            ylim([-1-sl 1+sl])
            xlabel('X location [m]')
            ylabel('Y location [m]')
            title('Inhomogenous response')
            axis square
            
        %Save the current frame into the video
        writeVideo(vidfile, getframe(gcf));
    end
end

close(vidfile);


%{
Animation - Critically damped case
%define filename for video
video_filename = 'Criticallydamped'

%create video file
vidfile = VideoWriter(video_filename, 'MPEG-4');
vidfile.FrameRate = 30;
open(vidfile);

%Initialize frame
nv = ceil(length(t)/10);
t_out = zeros(nv, 1);
x_hom = zeros(nv, 1);
x_inh = zeros(nv, 1);
count = 0;
sl = 0.25;


figure(11)
set(gcf,'Position',[15 50 1350 775])

for i = 1:1:length(t)
    %Plot a frame every 10s
    if mod(i-1,10) == 0
        count = count + 1;
        t_out(count) = t(i);
        %Trial 6 corresponds to user-defined underdamped case
        x_hom(count) = x_RK4_free(7,i);
        x_inh(count) = x_RK4_force(7,i);
        
        %Define the axes in the clockwise direction
        y1 = [x_hom(count)+sl x_hom(count)-sl x_hom(count)-sl x_hom(count)+sl];
        y2 = [x_inh(count)+sl x_inh(count)-sl x_inh(count)-sl x_inh(count)+sl];
        x = [-sl -sl sl sl];
        
        subplot(1,2,1)
            fill(x, y1, 'b')
            xlim([-1-sl 1+sl])
            ylim([-1-sl 1+sl])
            xlabel('X location [m]')
            ylabel('Y location [m]')
            title('Homogenous response')
            axis square
        subplot(1,2,2)
            fill(x, y2, 'b')
            xlim([-1-sl 1+sl])
            ylim([-1-sl 1+sl])
            xlabel('X location [m]')
            ylabel('Y location [m]')
            title('Inhomogenous response')
            axis square
            
        %Save the current frame into the video
        writeVideo(vidfile, getframe(gcf));
    end
end

close(vidfile);
%}


%{
Animation - Overdamped case
%define filename for video
video_filename = 'Overdamped'

%create video file
vidfile = VideoWriter(video_filename, 'MPEG-4');
vidfile.FrameRate = 30;
open(vidfile);

%Initialize frame
nv = ceil(length(t)/10);
t_out = zeros(nv, 1);
x_hom = zeros(nv, 1);
x_inh = zeros(nv, 1);
count = 0;
sl = 0.25;


figure(12)
set(gcf,'Position',[15 50 1350 775])

for i = 1:1:length(t)
    %Plot a frame every 10s
    if mod(i-1,10) == 0
        count = count + 1;
        t_out(count) = t(i);
        %Trial 6 corresponds to user-defined underdamped case
        x_hom(count) = x_RK4_free(8,i);
        x_inh(count) = x_RK4_force(8,i);
        
        %Define the axes in the clockwise direction
        y1 = [x_hom(count)+sl x_hom(count)-sl x_hom(count)-sl x_hom(count)+sl];
        y2 = [x_inh(count)+sl x_inh(count)-sl x_inh(count)-sl x_inh(count)+sl];
        x = [-sl -sl sl sl];
        
        subplot(1,2,1)
            fill(x, y1, 'b')
            xlim([-1-sl 1+sl])
            ylim([-1-sl 1+sl])
            xlabel('X location [m]')
            ylabel('Y location [m]')
            title('Homogenous response')
            axis square
        subplot(1,2,2)
            fill(x, y2, 'b')
            xlim([-1-sl 1+sl])
            ylim([-1-sl 1+sl])
            xlabel('X location [m]')
            ylabel('Y location [m]')
            title('Inhomogenous response')
            axis square
            
        %Save the current frame into the video
        writeVideo(vidfile, getframe(gcf));
    end
end

close(vidfile);
%}




