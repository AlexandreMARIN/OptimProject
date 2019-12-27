%we test the function gradient() in the folder SQP
clearvars
addpath("../SQP")

x = [0;0];%we will compute the gradient at x
h = [0.01, 0.01];%steps
[fx, cx] = func_constraints(x);
[g, J] = gradient(x, fx, cx, h, @func_constraints);

%we display the gradient of the objective and
%the Jacobian matrix of the linear constaints
disp(g)
disp(J)