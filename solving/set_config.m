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
[theta, R, V, M, tspanout] = path_solver();
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
    [theta, R, V, M, tspanout] = path_solver();
    V_r = norm(V(end,:));
    deltaV = V_c - V_r;
    fprintf("iter : %d\ndeltaV : %f\nm_e : %s\ntheta : %s\n\n", iter, deltaV, join(string(m_e)), join(string(theta)))

    iter = iter + 1;
end


%specific times
t_ = cumsum(t_c);

close all
%we draw the Earth and the orbit to reach
angle = linspace(-pi/300, pi/6, 1000);
xEarth = R_t*cos(angle);
yEarth = R_t*sin(angle);
xOrbit = R_c*cos(angle);
yOrbit = R_c*sin(angle);
figure(1)
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

%we compute norms of velocity and position
nbpt = size(tspanout, 1);
Vmag = zeros(nbpt, 1);
Rmag = zeros(nbpt, 1);
for n = 1:nbpt
    Vmag(n) = norm(V(n, :));
    Rmag(n) = norm(R(n, :));
end

figure(2)
plot(tspanout, Rmag, [0, t_(3)], [R_c, R_c])
xlabel("time $t$", 'Interpreter', 'latex')
ylabel("distance from the origin")
title({"the distance from the center of the Earth", "with respect to the time"})
xticks(t_)
xticklabels(["t_1", "t_2", "t_3"])
legend("$\Vert R\Vert$", "$R_c$", 'Interpreter', 'latex', 'Location', 'southeast')

figure(3)
plot(tspanout, Vmag, [0, t_(3)], [V_c, V_c], 'r--')
xlabel("time $t$", 'Interpreter', 'latex')
ylabel("velocity")
title("evolution of the velocity with respect to the time")
xticks(t_)
xticklabels(["t_1", "t_2", "t_3"])
legend("$\Vert V\Vert$", "$V_c$", 'Interpreter', 'latex', 'Location', 'southeast')

figure(4)
plot(tspanout, M)
xlabel("time $t$", 'Interpreter', 'latex')
ylabel("mass $M$", 'Interpreter', 'latex')
title("evolution of the mass with respect to the time")
xticks(t_)
xticklabels(["t_1", "t_2", "t_3"])
yticks(M_i([4 3 2 1]))
yticklabels(["m_u", "M_{i3}", "M_{i2}", "M_{i1}"])

fprintf("\nResults:\n\nV_p : %f\nfinal velocity : %f\ntheta : %s\nm_e : %s\n", V_p, Vmag(end), join(string(theta)), join(string(m_e)))