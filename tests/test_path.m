%we compute the path of the rocket
addpath("../simulator")
clearvars

%parameters of the problem
global R_t M_i alpha t_c v_e m_e k_ m_u

R_t = 6378137;
v_e = [2647.2 ; 2922.4 ; 4344.3];
k_ = [0.1101 ; 0.1532 ; 0.2154];
m_u = 1700;
m_e = [145349 ; 31215 ; 7933];
alpha = [15 ; 10 ; 10];

%we compute M_i and t_c
M_i = [0 ; 0 ; 0 ; m_u];
t_c = [0 ; 0 ; 0];
for j=[3 2 1]
    M_i(j) = M_i(j+1) + (1+k_(j))*m_e(j);
    t_c(j) = (m_e(j)*v_e(j))/(alpha(j)*M_i(j));
end

theta = [pi/10 ; pi/20 ; pi/50 ; pi/10];

%we solve the ODE
[R, V, M, tspanout] = rocketpath(theta);

%we plot graphs

close all
%we draw the Earth
angle = linspace(-pi/300, pi/10, 1000);
xEarth = R_t*cos(angle);
yEarth = R_t*sin(angle);
plot(xEarth, yEarth)
hold on
axis equal

%we draw the path
plot(R(:,1), R(:,2))

title("Example : path")
xlabel("$x$", 'Interpreter', 'latex')
ylabel("$y$", 'Interpreter', 'latex')
legend("Earth", "rocket")

nbpt = size(tspanout, 1);
Vmag = zeros(nbpt, 1);
Rmag = zeros(nbpt, 1);
for n = 1:nbpt
    Vmag(n) = norm(V(n, :));
    Rmag(n) = norm(R(n, :));
end
figure(2)
plot(tspanout, Rmag)
xlabel("time $t$", 'Interpreter', 'latex')
ylabel("$\Vert R\Vert$", 'Interpreter', 'latex')
title({"the distance from the center of the Earth", "with respect to the time"})

figure(3)
plot(tspanout, Vmag)
xlabel("time $t$", 'Interpreter', 'latex')
ylabel("$\Vert V\Vert$", 'Interpreter', 'latex')
title("velocity with respect to the time")

figure(4)
plot(tspanout, M)
xlabel("time $t$", 'Interpreter', 'latex')
ylabel("mass $M$", 'Interpreter', 'latex')
title("mass with respect to the time")