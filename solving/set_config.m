clearvars
addpath("../simulator", "../SQP")
global R_t R_c V_p V_c M_i m_u t_c m_e v_e k_ alpha

%we prepare some values
data()

%see that assignment in the first flow chart
V_p = 1.2*V_c;

%tolerance
tol = 15;

maxnbiter = 5;%maximum number of iterations
iter = 1;

%iteration : stage problem -> configuration: masses, thrusts, durations
%  -> path problem -> real velocity -> deltaV -> V_p is updated
stage_solver();
for j=[3 2 1]
    M_i(j) = M_i(j+1) + (1+k_(j))*m_e(j);
    t_c(j) = (m_e(j)*v_e(j))/(alpha(j)*M_i(j));
end
[theta, R, V, ~] = path_solver();
V_r = norm(V(end,:));
deltaV = V_c - V_r;

fprintf("iter : 0\ndeltaV : %f\nm_e : %s\ntheta : %s\n\n", deltaV, join(string(m_e)), join(string(theta)))


while (iter <= maxnbiter) && (abs(deltaV) > tol)
    V_p = V_p + deltaV;
    stage_solver();
    for j=[3 2 1]
        M_i(j) = M_i(j+1) + (1+k_(j))*m_e(j);
        t_c(j) = (m_e(j)*v_e(j))/(alpha(j)*M_i(j));
    end
    [theta, R, V, ~] = path_solver();
    V_r = norm(V(end,:));
    deltaV = V_c - V_r;fprintf("iter : %d\ndeltaV : %f\nm_e : %s\ntheta : %s\n\n", iter, deltaV, join(string(m_e)), join(string(theta)))
    iter = iter + 1;
end


close all
%we draw the Earth and the orbit to reach
angle = linspace(-pi/300, pi/6, 1000);
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
title("The path for $\theta = $["+join(string(theta))+"]", 'Interpreter', 'latex')
xlabel("$x$", 'Interpreter', 'latex')
ylabel("$y$", 'Interpreter', 'latex')
legend("Earth", "orbit", "rocket's path", "final velocity")

fprintf("\nResults:\n\nV_p : %f\nfinal velocity : %f\ntheta : %s\nm_e : %s\n", V_p, norm(V(end,:)), join(string(theta)), join(string(m_e)))