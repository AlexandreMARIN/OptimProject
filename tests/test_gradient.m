addpath("../")

x = [0;0];
h = 0.01;

[g, J] = gradient(x, func(x), constraints(x), h, @func, @constraints);
disp(g)
disp(J)