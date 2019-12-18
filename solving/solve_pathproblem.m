%we search for a value for initialising theta
%we assume the first problem is solved
%so m_e, M_i are correctly initialized
addpath("../simulator", "../SQP")
clearvars

global R_t R_c V_c M_i alpha m_u k_ m_e t_c v_e
%we prepare some values
data();
%we load m_e
load("m_e.mat")

%we compute M_i and t_c
for j=[3 2 1]
    M_i(j) = M_i(j+1) + (1+k_(j))*m_e(j);
    t_c(j) = (m_e(j)*v_e(j))/(alpha(j)*M_i(j));
end

%init_theta = [0.35 ; 0.18 ; 0.12 ; 0.08];%->you can try
%initial value for theta
init_theta = [0.05 ; 0.01 ; 0.05 ; 0.01];
[R, V, M, tspanout] = rocketpath(init_theta);

close all
%we draw the Earth and the orbit to reach
angle = linspace(-pi/300, pi/10, 1000);
xEarth = R_t*cos(angle);
yEarth = R_t*sin(angle);
xOrbit = R_c*cos(angle);
yOrbit = R_c*sin(angle);
plot(xEarth, yEarth, 'r--', xOrbit, yOrbit, 'g--')
hold on
axis equal

%we draw the path and the final velocity
plot(R(:,1), R(:,2))
quiver(R(end,1), R(end,2), V(end,1), V(end,2), 100, "+r-", "ShowArrowHead", "on", "MarkerSize", 10)
title("Trial 1 : the path for $\theta = $["+join(string(init_theta))+"]", 'Interpreter', 'latex')
xlabel("$x$", 'Interpreter', 'latex')
ylabel("$y$", 'Interpreter', 'latex')
legend("Earth", "orbit", "rocket's path", "final velocity")

%we compute norms of velocity and position
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
title({"Trial 1 : the distance from the center of the Earth", "with respect to the time"})

figure(3)
plot(tspanout, Vmag)
xlabel("time $t$", 'Interpreter', 'latex')
ylabel("$\Vert V\Vert$", 'Interpreter', 'latex')
title("Trial 1 : evolution of the velocity with respect to the time")

figure(4)
plot(tspanout, M)
xlabel("time $t$", 'Interpreter', 'latex')
ylabel("mass $M$", 'Interpreter', 'latex')
title("Trial 1 : evolution of the mass with respect to the time")


tol = 0.01;
theta_inf = init_theta-0.1;
theta_sup = init_theta+0.1;
maxnbiter = 20;
maxnbcall = 400;
mindtheta = 1e-5;
mindf = 1e-4;
%here we optimise with respect to the two last components of theta
[theta34, ~, iter, normGradLag] = SQP(init_theta(3:4), @(t34)pathproblem([init_theta(1:2);t34]), tol, theta_inf(3:4), theta_sup(3:4), maxnbiter, maxnbcall, mindtheta, mindf);

init_theta = [init_theta(1:2) ; theta34];
[R, V, M, tspanout] = rocketpath(init_theta);

figure(5)
%we draw the Earth and the orbit to reach
plot(xEarth, yEarth, 'r--', xOrbit, yOrbit, 'g--')
hold on
axis equal
%we draw the path
plot(R(:,1), R(:,2))
quiver(R(end,1), R(end,2), V(end,1), V(end,2), 100, "+r-", "ShowArrowHead", "on", "MarkerSize", 10)
title("Trial 2 : the path for $\theta = $["+join(string(init_theta))+"]", 'Interpreter', 'latex')
xlabel("$x$", 'Interpreter', 'latex')
ylabel("$y$", 'Interpreter', 'latex')
legend("Earth", "orbit", "rocket's path", "final velocity")

tol = 1e-6;
mindtheta = 1e-9;
mindf = 1e-6;

%here we optimise with respect to the last three components of theta
[theta234, ~, iter, normGradLag] = SQP(init_theta(2:4), @(t234)pathproblem([init_theta(1);t234]), tol, init_theta(2:4)-0.5, init_theta(2:4)+0.5, maxnbiter, maxnbcall, mindtheta, mindf);

init_theta = [init_theta(1) ; theta234];
[R, V, M, tspanout] = rocketpath(init_theta);

figure(6)
%we draw the Earth and the orbit to reach
plot(xEarth, yEarth, 'r--', xOrbit, yOrbit, 'g--')
hold on
axis equal
%we draw the path
plot(R(:,1), R(:,2))
quiver(R(end,1), R(end,2), V(end,1), V(end,2), 100, "+r-", "ShowArrowHead", "on", "MarkerSize", 10)
title("Trial 3 : the path for $\theta = $["+join(string(init_theta))+"]", 'Interpreter', 'latex')
xlabel("$x$", 'Interpreter', 'latex')
ylabel("$y$", 'Interpreter', 'latex')
legend("Earth", "orbit", "rocket's path", "final velocity")

mindf = 1e-11;
mindtheta = 1e-11;
[theta, ~, iter, normGradLag] = SQP(init_theta, @pathproblem, tol, 0.99998*init_theta, 1.00002*init_theta, maxnbiter, maxnbcall, mindtheta, mindf);
[R, V, M, tspanout] = rocketpath(theta);

figure(7)
%we draw the Earth and the orbit to reach
plot(xEarth, yEarth, 'r--', xOrbit, yOrbit, 'g--')
hold on
axis equal
%we draw the path
plot(R(:,1), R(:,2))
quiver(R(end,1), R(end,2), V(end,1), V(end,2), 100, "+r-", "ShowArrowHead", "on", "MarkerSize", 10)
title("Trial 4 : the path for $\theta = $["+join(string(theta))+"]", 'Interpreter', 'latex')
xlabel("$x$", 'Interpreter', 'latex')
ylabel("$y$", 'Interpreter', 'latex')
legend("Earth", "orbit", "rocket's path", "final velocity")

figure(8)
plot(tspanout, Vmag, [tspanout(1), tspanout(end)], [V_c V_c])
legend("velocity", "$V_c$", 'Interpreter', 'latex', 'Location', 'southeast')
xlabel("time $t$", 'Interpreter', 'latex')
ylabel("$\Vert V\Vert$", 'Interpreter', 'latex')
title("Trial 4 : evolution of the velocity with respect to the time")