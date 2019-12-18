%this script solves the stage problem
%You must initialize correctly V_p
addpath("../SQP")
clearvars

%parameters of the problem
global v_e m_u k_ V_p

%we prepare some values given by the problem
data();
V_p = V_p*1.2;%good

%lower and upper bounds, initial value given by the analytical solution
m_inf = [25000 ; 12000 ; 5000];
m_sup = [28000 ; 14000 ; 6000];
m0 = 1e4*[2.6981...
         ;1.2704...
         ;0.5516];

mindm = 1;
mindf = 0.1;
maxnbcall = 300;
maxnbiter = 20;
tol = 0.01;

[m_e, lambda, iter, normGradLag, f_m, c_m] = SQP(m0, @stageproblem, tol, m_inf, m_sup, maxnbiter, maxnbcall, mindm, mindf);

fprintf("f(m) = %f\n\nc(m) = %f\n\n", f_m, c_m)

%we save the value of m_e to reuse it
save("m_e.mat", "m_e")