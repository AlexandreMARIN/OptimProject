function dydt = odefunc(t, y, alpha, M_i, v_e, theta)
%y := R, V, M
global R_t

mu = 3.986e14;
c_x = 0.1;
rho0 = 1.225;
H = 7e3;

rho = rho0 * exp( (R_t - norm(y(1:2)))/H );
Tmag = alpha*M_i;

R = norm(y(1:2));
V = norm(y(3:4));

dydt = zeros(5, 1);

dydt(1) = y(3);
dydt(2) = y(4);


W = -(mu*y(5)/norm(y(1:2))^3)*y(1:2);
D = -c_x*rho*norm(y(3:4))*y(3:4);

e_r = y(1:2)/norm(y(1:2));
e_h = [-e_r(2) ; e_r(1)];
gamma = asin(y(1:2)'*y(3:4)/(R*V));
angle = theta + gamma;
u = cos(angle)*e_h + sin(angle)*e_r;
T = Tmag*u;

dydt(3:4) = (T + W + D)/y(5);

dydt(5) = -Tmag/v_e;

end

