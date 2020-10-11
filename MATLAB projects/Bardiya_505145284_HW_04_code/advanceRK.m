function [ y ] = advanceRK (y, dt, method)
%The following function performs the first, second or fourth order
%Runge-Kutta method according to user-input via advancing the discretized
%solution by 1 time-step. The current amount of y is used to calculate the
%value at the next time-step and is finally returned as output
t_half = 2.45;
switch(method)
    %First Runge-Kutta Method
    case (1)
        c1 = dt*(-log(2)/t_half)*y;
        y_next = y + c1;
    %Second Runge-Kutta Method 
    case (2)
        c1 = dt*(-log(2)/t_half)*y;
        c2 = dt*(-log(2)/t_half)*(y + 0.5*c1);
        y_next = y + c2;
    %Fourth Runge-Kutta Method
    case (4)
        c1 = dt*(-log(2)/t_half)*y;
        c2 = dt*(-log(2)/t_half)*(y + 0.5*c1);
        c3 = dt*(-log(2)/t_half)*(y + 0.5*c2);
        c4 = dt*(-log(2)/t_half)*(y + c3);
        y_next = y + c1/6 + c2/3 + c3/3 + c4/6;
    otherwise
        fprintf('Error: You must choose 1,2 or 4 to select a particular Runge-Kutta method.\n');
end
%Update value of y
y = y_next;



