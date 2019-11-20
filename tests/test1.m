addpath("../")


x0 = [-1 ; 2 ; 1 ; -2 ; -2];

x_star = [-1.2366 ; 2.4616 ; 1.1911 ; -0.2144 ; -1.6165];

f_star = 28.4974;

fprintf("f(x_star) = %f\n\n", f_test1(x_star))
fprintf("c(x_star) = %f\n\n", c_test1(x_star))
x_inf = x_star-0.3;
x_sup = x_star+0.3;

[x,lambda,k] = SQP(x0, @f_test1, @c_test1, 0.001, 10000,x_inf,x_sup);
