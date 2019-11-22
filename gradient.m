function [g, J] = gradient(x, fx, cx, h, f, c)
%Input:
%x : point of size n*1
%fx, cx : values of f and c at x
%h : step
%f : real-valued function of n variables
%c : function which gives the m constraints
%Output :
%g : gradient of f at x
%J : the Jacobian matrix of c at x

[n, ~] = size(x);

[m, ~] = size(cx);

g = zeros(n, 1);
J = zeros(m, n);

for i=1:n
    xh = [x(1:i-1) ; x(i) + h ; x(i+1:n)];
    g(i) = (f(xh) - fx)/h;%df/dx_i
    J(:, i) = (c(xh) - cx)/h;
end


end