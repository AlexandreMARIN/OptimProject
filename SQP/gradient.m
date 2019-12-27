function [g, J] = gradient(x, fx, cx, h, problem)
%Input:
%x : point of size n*1
%fx, cx : values of f and c at x
%h : steps, i.e. a vector giving a step for each component of x
%problem : it must be the same thing as in the SQP function
%Output :
%g : gradient of f at x
%J : the Jacobian matrix of c at x

n = size(x, 1);
m = numel(problem{2});%number of constraints

g = zeros(n, 1);

if m == 0
    J = zeros(0, n);
    for i=1:n
        xh = [x(1:i-1) ; x(i) + h(i) ; x(i+1:n)];
        [fxh, cxh] = problem{1}(xh);
        cxh = cxh(problem{2});
        g(i) = (fxh - fx)/h(i);%df/dx_i
    end
    return
end

J = zeros(m, n);

for i=1:n
    xh = [x(1:i-1) ; x(i) + h(i) ; x(i+1:n)];
    [fxh, cxh] = problem{1}(xh);
    cxh = cxh(problem{2});
    g(i) = (fxh - fx)/h(i);%df/dx_i
    J(:, i) = (cxh - cx)/h(i);
end


end