function [H2] = compute_H(H1, x1, fx1, cx1, x2, fx2, cx2, lambda, f, c)
%returns a symmetric positive definite n-by-n matrix H2
%Input:
%H1 : a n-by-n symmetric positive definite matrix
%x1 : the former value of x
%fx1, cx1 : values of f and c at x1
%x2 : the current value of the estimate of a minimizer x
%fx2, cx2 : values of f and c at x2
%lambda : the current multipliers of Lagrange
%f : the objective function
%c : the constraints
d = x2 - x1;

h = 0.01;%step for the gradient

[g2, J2] = gradient(x2, fx2, cx2, h, f, c);

y = g2 + J2' * lambda;

[g, J] = gradient(x1, fx1, cx1, h, f, c);

y = y - (g + J' * lambda);


H2 = BFGS(H1, y, d);

end