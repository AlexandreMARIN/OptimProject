addpath("../")
clearvars

x0 = [-1 ; 2 ; 1 ; -2 ; -2];

x_star = [-1.2366 ; 2.4616 ; 1.1911 ; -0.2144 ; -1.6165];

f_star = 28.4974;

fprintf("f(x_star) = %f\n\n", f_test1(x_star))
fprintf("c(x_star) = %f\n\n", c_test1(x_star))
x_inf = x_star - 2;
x_sup = x_star + 2;
tol = 1e-7;
maxnbiter = 50;
maxnbcall = 100;
mindx = 1e-8;
mindf = 1e-10;
[x, lambda, iter, normGradLag] = SQP(x0, @f_test1, @c_test1, tol, x_inf, x_sup, maxnbiter, maxnbcall, mindx, mindf);
