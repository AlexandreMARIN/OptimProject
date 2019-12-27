%test case 1
addpath("../SQP")
clearvars

%the starting point
x0 = [-1 ; 2 ; 1 ; -2 ; -2];

%the searched minimiser
x_star = [-1.2366 ; 2.4616 ; 1.1911 ; -0.2144 ; -1.6165];

f_star = 28.4974;%value of the objective at
[fxstar, cxstar] = fc_test1(x_star);
%we check the minimiser
fprintf("f(x_star) = %f\n\n, c(x_star) = %f\n\n", fxstar, cxstar)

%parameters for the SQP function
x_inf = x_star - 0.3;
x_sup = x_star + 0.3;
tol = 1e-7;
maxnbiter = 10;
maxnbcall = 60;
mindx = 1e-8;
mindf = 1e-10;

%we optimize
[x, lambda, iter, normGradLag, fx, cx, rho, nbcall] = SQP(x0, @fc_test1, tol, x_inf, x_sup, maxnbiter, maxnbcall, mindx, mindf, [true true]);
