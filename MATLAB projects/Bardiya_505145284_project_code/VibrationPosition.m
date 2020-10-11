function [ x ] = VibrationPosition(x0,m,k,c,f,dt,type)
%VibrationPosition function
%x0 = current conditions (position, velocity)
%m = mass
%k = spring constant
%c = damping constant
%f = force function which is 1*3 vector where the force is determined
%according to the particular Runge-Kutta method
%dt = time-step
%type = either Runge-Kutta 1/2/4 for the numerical integration method

x_k = x0(1);
v_k = x0(2);

%Natural Frequency
w_n = sqrt(k/m);
%Damping ratio
xi = c/(2*m*w_n);

switch (type)
    %Forward Euler Method (Runge-Kutta 1)
    case (1)
        %Velocity - first derivative of position
        dxdt = v_k;
        %accelaration - second derivative of position
        dvdt = -2*xi*w_n*v_k - w_n^2*x_k + f(1);
        
        %Calculate the values for the next iteration
        x_kp1 = x_k + dt*dxdt;
        v_kp1 = v_k + dt*dvdt;
        
        %Update the initial conditions
        x = [x_kp1 v_kp1];
        
    %Runge-Kutta 2
    case (2)
        %Velocity - first derivative of position
        dxdt = v_k;
        %accelaration - second derivative of position
        dvdt = -2*xi*w_n*v_k - w_n^2*x_k + f(1);
        
        %Implementing 2nd order Runge-Kutta
        c_x1 = dt*dxdt;
        c_v1 = dt*dvdt;
        c_x2 = dt*(dxdt + 0.5*c_v1);
        c_v2 = dt*(-2*xi*w_n*(v_k + 0.5*c_v1) - w_n^2*(x_k + 0.5*c_x1) + f(2));
        
        %Calculate the values for the next iteration
        x_kp1 = x_k + c_x2;
        v_kp1 = v_k + c_v2;
        
        %Update the initial conditions
        x = [x_kp1 v_kp1];
        
    %Runge-Kutta 4
    case (4)
        %Velocity - first derivative of position
        dxdt = v_k;
        %accelaration - second derivative of position
        dvdt = -2*xi*w_n*v_k - w_n^2*x_k + f(1);
        
        %Implementing 4th order Runge-Kutta
        c_x1 = dt*dxdt;
        c_v1 = dt*dvdt;
        c_x2 = dt*(dxdt + 0.5*c_v1);
        c_v2 = dt*(-2*xi*w_n*(v_k + 0.5*c_v1) - w_n^2*(x_k + 0.5*c_x1) + f(2));
        c_x3 = dt*(dxdt + 0.5*c_v2);
        c_v3 = dt*(-2*xi*w_n*(v_k + 0.5*c_v2) - w_n^2*(x_k + 0.5*c_x2) + f(2));
        c_x4 = dt*(dxdt + c_v3);
        c_v4 = dt*(-2*xi*w_n*(v_k + c_v3) - w_n^2*(x_k + c_x3) + f(3));
        
        %Calculate the values for the next iteration
        x_kp1 = x_k + (c_x1 + 2*c_x2 + 2*c_x3 + c_x4)/6;
        v_kp1 = v_k + (c_v1 + 2*c_v2 + 2*c_v3 + c_v4)/6;
        
        %Update the initial conditions
        x = [x_kp1 v_kp1];
        
    otherwise
        fprintf('Error: Please select 1,2 or 4 to choose a Runge-Kutta method.\n');
end
