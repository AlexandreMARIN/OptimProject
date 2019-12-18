clearvars
addpath("../SQP")

x = [0;0];
h = [0.01, 0.01];
[fx, cx] = func_constraints(x);
[g, J] = gradient(x, fx, cx, h, @func_constraints);
disp(g)
disp(J)