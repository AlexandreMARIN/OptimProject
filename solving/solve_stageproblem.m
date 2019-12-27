%this script solves the stage problem
%You must initialize correctly V_p
addpath("../SQP")
clearvars

%parameters of the problem
global V_p

%we prepare some values given by the problem
data();
V_p = V_p*1.2;%good

%lower and upper bounds, initial value given by the analytical solution
m_inf = [25000 ; 12000 ; 5000];
m_sup = [28000 ; 14000 ; 6000];

%the solution of the first problem with V_p = 1.2*V_c is
%m0 = 1e4*[2.6981...
%         ;1.2704...
%         ;0.5516];
m0 = 1e4*[2.7...
         ;1.3...
         ;0.6];

%parameters for the function SQP()
mindm = 0.01;
mindf = 0.001;
maxnbcall = 180;
maxnbiter = 20;
tol = 0.001;

%we optimize
[m_e, lambda, iter, normGradLag, f_m, c_m] = SQP(m0, @stageproblem,...
    tol, m_inf, m_sup, maxnbiter, maxnbcall, mindm, mindf);

%we display the values of the objective and of the constraint
fprintf("m_e: %s\n\nf(m_e) = %f\n\nc(m_e) = %f\n\n",...
        join(string(m_e), ", "), f_m, c_m)

%we save the value of m_e to reuse it
save("m_e.mat", "m_e")