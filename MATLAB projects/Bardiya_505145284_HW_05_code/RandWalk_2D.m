function [x, y] = RandWalk_2D(x0, y0, BC)
%A random-walk function that takes the initial position of the particle
%(x0, y0) and computes the directional movement of the particle; each
%possible direction having an equal probablity of occuring. Takes into
%account boundary conditions specified in vector BC

%Create a random number for determining the directiong of movement
r = rand;

%Move North
if r <= 0.2
    x = x0;
    y = y0 + 1;
    %check boundary condition
    if y >= BC(1)
        y = BC(1);
    end 
%Move South
elseif r > 0.2 && r <= 0.4
    x = x0;
    y = y0 - 1;
    %check boundary condition
    if y <= BC(2)
        y = BC(2);
    end
%Move West
elseif r > 0.4 && r <= 0.6
    x = x0 - 1;
    y = y0;
    %check boundary condition
    if x <= BC(3)
        x = BC(3);
    end
%Move East
elseif r > 0.6 && r <= 0.8
    x = x0 + 1;
    y = y0;
    %check boundary condition
    if x >= BC(4)
        x = BC(4);
    end
%No Movement
elseif r > 0.8
    x = x0;
    y = y0;
end  
    