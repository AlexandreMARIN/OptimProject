addpath("../")
global R_t

R_t = 6378137;
R0 = [R_t ; 0];
theta0 = -pi/100;
V0mag = 100;
V0 = V0mag*[cos(theta0) ; sin(theta0)];
M0 = 208691;

alpha = 15;
v_e = 2600;
m_e = 145349;
tf = m_e*v_e/(alpha*M0);
t0 = 0;


%for the Cauchy problem
y0 = [R0 ; V0 ; M0];
tspan = linspace(t0, tf, 1000);

%we solve the ODE
[t, y] = ode45(@(t, y) odefunc(t, y, alpha, M0, v_e, theta0), tspan, y0);

close all
%we draw the Earth
angle = linspace(-pi/200, pi/200, 1000);
xEarth = R_t*cos(angle);
yEarth = R_t*sin(angle);
plot(xEarth, yEarth)
hold on
axis equal

%we draw the path
plot(y(:,1), y(:,2))

title("Example: the first part of the path")
xlabel("x")
ylabel("y")
legend("Earth", "rocket")
