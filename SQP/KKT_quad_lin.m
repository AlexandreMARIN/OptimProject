function [x, y] = KKT_quad_lin(A, b, C, d, getWarn)
%Input:
%A is a n-by-n symmetric positive definite matrix
%b is a vector of size n
%C is a m-by-n matrix of rank m
%d is a vector of size m
%getWarn : a Boolean value : true enables this function to raise
%a warning if necessary
%Output:
%that function minimizes the function
%0.5*<Ax, x> + <b, x>
%with the constraint
%C*x = d
%x : the minimizer
%y : the multipliers of Lagrange
[m, n] = size(C);
M = [A C' ; C zeros(m)];
opts.SYM = true;

if getWarn
    z = linsolve(M, [-b ; d], opts);
else
    [z, ~] = linsolve(M, [-b ; d], opts);
end

x = z(1:n);
if m == 0
    y = [];
else
    y = z(n+1:n+m);
end

end