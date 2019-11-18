function [x, y] = KKT_quad_lin(A, b, C, d)
%Input:
%A is a n-by-n symmetric positive definite matrix
%b is a vector of size n
%C is a m-by-n matrix of rank m
%d is a vector of size m
%Output:
%that function minimizes the function
%0.5*<Ax, x> + <b, x>
%with the constraint
%C*x = d
%x : the minimizer
%y : the multipliers of Lagrange
[m, n] = size(C);
M = [A C' ; C zeros(m)];

z = linsolve(M, [-b ; d]);

x = z(1:n);
y = z(n+1:m);

end