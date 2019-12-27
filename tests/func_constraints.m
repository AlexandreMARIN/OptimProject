function [fx, cx] = func_constraints(x)
%use this function with the file 'test_gradient.m'
%x must be a 2D real vector

fx = cos(x(1)) + exp(x(2));

A = [1 0 ; 1 1 ; 0 1];

b = [0 ; 1 ; 0];

cx = A*x + b;

end