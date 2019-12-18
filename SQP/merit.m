function [F] = merit(fx, cx, rho)
%Input:
%fx, cx : evaluation of f and c at a point x
%rho : a real number
%Output:
%F : the evaluation of the merit function at x

F = fx + rho*norm(cx, 1);


end