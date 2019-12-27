%we search for a good initial value of theta
%we assume the first problem is solved
%so m_e, M_i are or will be correctly initialized

addpath("../simulator", "../SQP")
clearvars
close all
global R_t R_c V_c M_i alpha k_ m_e t_c v_e

tolV = 50;%tolerance about |V_c - final velocity|

%we set some parameters for SQP()
tol = 0.5;
maxnbiter = 5;
maxnbcall = 50;
mindtheta = 1e-4;
mindf = 1e-3;

%intervals for each component of theta
interv1 = linspace(0, pi/5, 5);
interv2 = linspace(0, pi/9, 5);
interv3 = linspace(0, pi/9, 5);
interv4 = linspace(0, pi/9, 5);

%we prepare some values
data()
%we load m_e
load("m_e.mat")

%we compute M_i and t_c
for j=[3 2 1]
    M_i(j) = M_i(j+1) + (1+k_(j))*m_e(j);
    t_c(j) = (m_e(j)*v_e(j))/(alpha(j)*M_i(j));
end

%we draw the Earth and the orbit to reach
figure(1)
angle = linspace(-pi/300, pi/10, 1000);
xEarth = R_t*cos(angle);
yEarth = R_t*sin(angle);
xOrbit = R_c*cos(angle);
yOrbit = R_c*sin(angle);
plot(xEarth, yEarth, 'r--', xOrbit, yOrbit, 'g--')
hold on
axis equal

leg = string(zeros(1, 1250));%the legend of the graph
theta_nb = 0;%number of accepted values for theta

for theta1 = interv1
    for theta2 = interv2
        for theta3 = interv3
            for theta4 = interv4
                %first, we optimize according to the constraint about
                %altitude and according to the first two components of
                %theta
                theta12_inf = [theta1 ; theta2] - 0.005;
                theta12_sup = [theta1 ; theta2] + 0.005;
                theta12 = SQP([theta1 ; theta2], {@(t12) pathproblem([t12 ; theta3 ; theta4]) ; 1}, tol, theta12_inf, theta12_sup, maxnbiter, maxnbcall, mindtheta, mindf);
                if any(isnan(theta12))
                    continue
                end
                
                %then, we optimize according to the constraint about
                %the velocity and according to the last two components of
                %theta
                theta34_inf = [theta3 ; theta4] - 0.1;
                theta34_sup = [theta3 ; theta4] + 0.1;
                theta34 = SQP([theta3 ; theta4], {@(t34) pathproblem([theta12 ; t34]) ; 2}, tol, theta34_inf, theta34_sup, maxnbiter, maxnbcall, mindtheta, mindf);
                if any(isnan(theta34))
                    continue
                end
                
                %at last, we optimize normally
                theta = [theta12 ; theta34];
                theta_inf = [theta(1:2)-0.002 ; theta(3:4)-0.1];
                theta_sup = [theta(1:2)+0.002 ; theta(3:4)+0.1];
                [theta, ~, ~, normGradLag] = SQP(theta, @pathproblem, tol, theta_inf, theta_sup, maxnbiter, maxnbcall, mindtheta, mindf);
                if any(isnan(theta))
                    continue
                end
                
                %we select only the good values
                [R, V] = rocketpath(theta);
                deltaV = abs(V_c - norm(V(end,:)));
                if ( normGradLag > tol) || (deltaV > tolV)
                    continue
                end
                fprintf("theta: %s\ndeltaV : %f\n\n", join(string(theta), ", "), deltaV)
                
                %we draw the path and the final velocity
                plot(R(:,1), R(:,2))
                theta_nb = theta_nb + 1;
                leg(2*theta_nb-1) = sprintf("Path %d", theta_nb);
                leg(2*theta_nb) = sprintf("Final velocity %d", theta_nb);
                quiver(R(end,1), R(end,2), V(end,1), V(end,2), 100, "+r-", "ShowArrowHead", "on", "MarkerSize", 10)
            end
        end
    end
end

leg = ["Earth", "orbit", leg(1:2*theta_nb)];
legend(leg, 'Location', 'southeast')
xlabel("$x$", 'Interpreter', 'latex')
ylabel("$y$", 'Interpreter', 'latex')
title({"Some trials:", "which initial value of theta is good ?"})
