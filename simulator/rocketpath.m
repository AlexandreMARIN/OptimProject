function [R, V, M, tspanout] = rocketpath(theta)
%Input:
%theta : a four-value vector
%Output:
%R : position of the rocket
%V : the velocity of the rocket
%M : the mass of the rocket
%tspanout : times where R, V and M are evaluated


%data from the problem
global R_t
global M_i
global alpha
global v_e
global t_c

%results : [R, V, M]
res = [];

%for the Cauchy problem
R0 = [R_t ; 0];%initial position
V0mag = 100;
V0 = V0mag*[cos(theta(1)) ; sin(theta(1))];%initial velocity
sol = [R0 ; V0 ; 0]';%initial values of the solution of the ODE
tspan = [0, 0];%we integrate onto this time interval
tspanout = [];

for j = 1:3
    %we update the initial value of the derivative
    %and we change the time interval
    tspan(1) = tspan(2);
    tspan(2) = tspan(2) + t_c(j);
    y0 = [sol(end,1:4)' ; M_i(j)];
    %we solve the ODE
    [t, sol] = ode45(@(t, y) odefunc(t, y, alpha(j), M_i(j), v_e(j), theta(j+1)), tspan, y0);
    res = [res(:, :) ; sol];
    tspanout = [tspanout(:); t];
end

R = res(:, 1:2);
V = res(:, 3:4);
M = res(:, 5);

end

