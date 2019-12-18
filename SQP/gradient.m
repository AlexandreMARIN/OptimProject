function [g, J] = gradient(x, fx, cx, h, problem)
%Input:
%x : point of size n*1
%fx, cx : values of f and c at x
%h : steps, i.e. a vector giving a step for each component of x
%problem : a function which returns something of the form [f(x), c(x)]
%where
%  f is a real-valued function of n variables
%  c is a function which gives the m constraints
%Output :
%g : gradient of f at x
%J : the Jacobian matrix of c at x

[n, ~] = size(x);

[m, ~] = size(cx);

g = zeros(n, 1);
J = zeros(m, n);

for i=1:n
    xh = [x(1:i-1) ; x(i) + h(i) ; x(i+1:n)];
    [fxh, cxh] = problem(xh);
    g(i) = (fxh - fx)/h(i);%df/dx_i
    J(:, i) = (cxh - cx)/h(i);
end


end