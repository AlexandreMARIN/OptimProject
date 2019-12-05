function dydt = odefunc(t, y, alpha, M_i, v_e, theta)
%this function must be used with ode45()
%Input:
%t, y : for ode45
%here, y represents [R, V, M]
%other parameters : given by the problem

global R_t

mu = 3.986e14;
c_x = 0.1;
rho0 = 1.225;
H = 7e3;

%distance between the rocket and the center of the Earth
R = norm(y(1:2));
%velocity of the rocket
V = norm(y(3:4));

rho = rho0 * exp( (R_t - R)/H );
Tmag = alpha*M_i;

dydt = zeros(5, 1);

%\dot{R} = V
dydt(1:2) = y(3:4);

%forces
W = -(mu*y(5)/R^3)*y(1:2);
D = -c_x*rho*V*y(3:4);

e_r = y(1:2)/R;
e_h = [-e_r(2) ; e_r(1)];
gamma = asin(y(1:2)'*y(3:4)/(R*V));
angle = theta + gamma;
u = cos(angle)*e_h + sin(angle)*e_r;
T = Tmag*u;

%Fundamental principle of the dynamics
dydt(3:4) = (T + W + D)/y(5);

%derivative of the mass
dydt(5) = -Tmag/v_e;

end

