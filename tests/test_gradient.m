addpath("../")

x = [0;0];
h = 0.01;

[g, J] = gradient(x, h, @func, @constraints);
disp(g)
disp(J)