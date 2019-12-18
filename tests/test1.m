addpath("../SQP")
clearvars

x0 = [-1 ; 2 ; 1 ; -2 ; -2];

x_star = [-1.2366 ; 2.4616 ; 1.1911 ; -0.2144 ; -1.6165];

f_star = 28.4974;
[fxstar, cxstar] = fc_test1(x_star);
fprintf("f(x_star) = %f\n\n", fxstar)
fprintf("c(x_star) = %f\n\n", cxstar)
x_inf = x_star - 0.3;
x_sup = x_star + 0.3;
tol = 1e-7;
maxnbiter = 50;
maxnbcall = 100;
mindx = 1e-8;
mindf = 1e-10;
[x, lambda, iter, normGradLag] = SQP(x0, @fc_test1, tol, x_inf, x_sup, maxnbiter, maxnbcall, mindx, mindf);
