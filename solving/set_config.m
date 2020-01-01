%here we search for an optimum configuration of the rocket
clearvars
addpath("../simulator", "../SQP")
global R_t H_c R_c V_p V_c M_i t_c m_e v_e k_ alpha

%we prepare some values
data()

%see that assignment in the first flow chart
V_p = 1.2*V_c;

%tolerance
tol = 10;

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

fprintf("iteration : %d\ndeltaV : %f\nm_e : %s\nV_p : %f\nfinal velocity : %f\ntotal mass of the rocket : %f\ntheta : %s\n\n",...
    iter, deltaV, join(string(m_e)), V_p, V_r, M_i(1), join(string(theta)))

iter = iter + 1;

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
    fprintf("iteration : %d\ndeltaV : %f\nm_e : %s\nV_p : %f\nfinal velocity : %f\ntotal mass of the rocket : %f\ntheta : %s\n\n",...
            iter, deltaV, join(string(m_e)), V_p, V_r, M_i(1), join(string(theta)))

    iter = iter + 1;
end

nbiter = iter - 1;%number of iterations to get the results

%%%%%%%%%%%%%%%%%%%%%%%%
%we present some results
%%%%%%%%%%%%%%%%%%%%%%%%

%specific times : the combustion of the stage j ends at t_(j)
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
quiver(R(end,1), R(end,2), V(end,1), V(end,2), 50, "+r-", "ShowArrowHead", "on", "MarkerSize", 10, 'MaxHeadSize', 1)
title({"The path for", "$\theta = $[ "+join(string(theta), ", ")+" ]"}, 'Interpreter', 'latex')
xlabel("$x$", 'Interpreter', 'latex')
ylabel("$y$", 'Interpreter', 'latex')
legend("Earth", "orbit", "path of the rocket", "final velocity", 'FontSize', 12)

%we compute norms of the velocity and of the position
nbpt = size(tspanout, 1);
Vmag = zeros(nbpt, 1);
Rmag = zeros(nbpt, 1);
for n = 1:nbpt
    Vmag(n) = norm(V(n, :));
    Rmag(n) = norm(R(n, :));
end

%we plot the altitude
figure(2)
plot(tspanout, Rmag - R_t, [0, t_(3)], [H_c, H_c])
xlabel("time")
ylabel("altitude")
title({"the altitude of the rocket (in meters)", "with respect to the time"})
xticks(t_)
xticklabels(["t_1", "t_2", "t_3"])
legend("$\Vert R\Vert\ -\ R_{t}$", "$H_c$", 'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 15)

%we display the velocity
figure(3)
plot(tspanout, Vmag, [0, t_(3)], [V_c, V_c], 'r--')
xlabel("time")
ylabel("velocity")
title({"the velocity of the rocket (in m/s)", "with respect to the time"})
xticks(t_)
xticklabels(["t_1", "t_2", "t_3"])
legend("$\Vert V\Vert$", "$V_c$", 'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 15)

%we display the mass
figure(4)
M_f = M_i(1:3) - m_e;%final masses
plot(tspanout, M, "-k")
hold on
for j = [2 3]
    plot([t_(j)-t_c(j), t_(j)], [M_i(j), M_i(j)], '--m')
end
for j = [1 2]
    plot([t_(j)-t_c(j), t_(j)], [M_f(j), M_f(j)], '--b')
end
xlabel("time")
ylabel("mass")
title("the mass of the rocket with respect to the time")
xticks(t_)
xticklabels(["t_1", "t_2", "t_3"])
yticks(sort([M_i(1:3) ; M_f]))
yticklabels(["M_{f3}", "M_{i3}", "M_{f2}", "M_{i2}", "M_{f1}", "M_{i1}"])
ylim([M_f(3), M_i(1)])

%we display some information
fprintf("\nResults acquired in %d iterations:\n\nm_e : %s\nV_p : %f\nfinal velocity : %f\ntotal mass of the rocket : %f\ntheta : %s\n",...
    nbiter, join(string(m_e)), V_p, Vmag(end), M_i(1), join(string(theta)))