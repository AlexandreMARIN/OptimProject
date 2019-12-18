function [H2] = compute_H(H1, d, gxlag, prev_g, prev_A, lambda)
%returns a symmetric positive definite n-by-n matrix H2
%Input:
%H1 : a n-by-n symmetric positive definite matrix
%d : equals x2 - x1 where
%  x1 is the former estimate of the minimiser
%  x2 is the current estimate of the minimiser
%gxlag : value of the gradient with respect ro x, 
%of the Lagrangian of the objective, at x2
%prev_g : the previous value of the gradient of the objective (at x1)
%prev_A : the previous value of the Jacobian matrix of the constraints
%(at x1)
%lambda : the current multipliers of Lagrange

y = gxlag - (prev_g + prev_A'*lambda);

H2 = BFGS(H1, y, d);

end