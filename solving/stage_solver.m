function [] = stage_solver()
%it solves the stage problem
%Don't forget to initialize V_p correctly!
global m_e

%parameters for the function SQP()
mindm = 1;
mindf = 0.1;
maxnbcall = 80;
maxnbiter = 10;
tol = 0.01;

startnb = 50;%number of starting points we test

%intervals giving possible values of the components of m_e
%these intervals are described with their minimums (in the first column)
%and their lengths (in the second column)
%it is a neighbourhood of the solution of the problem with V_p=1.2*V_c
interv_m = [24000, 5000 ...
           ;10000, 4000 ...
           ;5000, 1000];

min_ngl = Inf;%minimum norm of the gradient of the Lagrangian function
%lower and upper bounds
m_inf = interv_m(:,1);
m_sup = interv_m(:,1)+interv_m(:,2);

%here we test several starting points and we select the estimate giving
%the least value for normGradLag
for i = 1:startnb
    init_m = interv_m(:,1) + rand(3, 1).*interv_m(:,2);
    [m, ~, ~, normGradLag] = SQP(init_m, @stageproblem, tol, m_inf, m_sup, maxnbiter, maxnbcall, mindm, mindf);
    if normGradLag < tol
        m_e = m;
        break
    end
    if normGradLag < min_ngl
        min_ngl = normGradLag;
        m_e = m;
    end
end

end