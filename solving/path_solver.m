function [theta, R, V, M, tspanout] = path_solver()
%This function is supposed to be used to get the optimum configuration
%We assume that the required data (like t_c) are correctly initialised
%Output:
%theta : an estimate of the solution of the path problem
%R, V, M, tspanout : position, velocity, mass of the rocket and times where
%R, V and M are evaluated

%we set some parameters for the SQP function
maxnbiter = 10;
maxnbcall = 100;
mindtheta = 1e-5;
mindf = 1e-4;
tol = 1e-4;

startnb = 50;%number of trials

%intervals giving possible values of the components of theta
%these intervals are described with their minimums (in the first column)
%and their lengths (in the second column)
%it is a neighbourhood of a solution of the initial problem with V_p=1.2V_c
interv_theta = [0.04, 0.013...
               ;-0.014, 0.009...
               ;0.26, 0.09...
               ;0.08, 0.08];

min_ngl = Inf;%minimum norm of the gradient of the Lagrangian function
%lower and upper bounds
theta_inf = interv_theta(:,1);
theta_sup = interv_theta(:,1)+interv_theta(:,2);
theta = [NaN ; NaN ; NaN ; NaN];

%here we test several starting points and we select the estimate giving
%the least value for normGradLag
for i = 1:startnb
    init_theta = interv_theta(:,1) + rand(4, 1).*interv_theta(:,2);
    [res, ~, ~, normGradLag] = SQP(init_theta, @pathproblem, tol, theta_inf, theta_sup, maxnbiter, maxnbcall, mindtheta, mindf);
    if normGradLag < tol
        theta = res;
        break
    end
    if normGradLag < min_ngl
        theta = res;
        min_ngl = normGradLag;
    end
end

%here we optimise with respect to the last two components of theta
theta34 = SQP(theta(3:4), @(t34)pathproblem([theta(1:2);t34]), tol, theta_inf(3:4), theta_inf(3:4), maxnbiter, maxnbcall, mindtheta, mindf);
theta = [theta(1:2) ; theta34];

%here we optimise with respect to the last three components of theta
theta234 = SQP(theta(2:4), @(t234)pathproblem([theta(1);t234]), tol, theta_inf(2:4), theta_sup(2:4), maxnbiter, maxnbcall, mindtheta, mindf);
theta = [theta(1) ; theta234];

[R, V, M, tspanout] = rocketpath(theta);

end