function [x, fx, cx, nbcall, dx] = Armijo(d, x, x0, fx0, gx0, cx0, c1, f, c, rho0, nbiter, nbcall)
%Input:
%d : d = x - x0, direction of descent
%x : the approximation of the minimiser given by the Newton's algorithm
%x0 : the former approximation of the minimiser
%fx0, gx0, cx0 : resp. values of f, gradient of f and c at x0
%c1 : for the Armijo's algorithm
%f, c : objective function and constraints
%rho0 : for the merit function, initial value of rho
%nbiter : number of iterations
%nbcall : number of calls for f, c
%
%Output:
%x : the new estimate of the minimiser
%fx, cx : f(x), c(x)
%nbcall : number of calls for f, c
%dx : = x - x0

%we compute the derivative of the merit function, with respect to the direction d
%we multiply it by c1
rho = rho0;
normcx0 = norm(cx0, 1);
gx0d = gx0'*d;
c1F_d = gx0d - rho*normcx0;
%we must have a direction of descent
%so we increase the value of rho
while c1F_d >= 0
    rho = rho * 10;
    c1F_d = gx0d - rho*normcx0;
end
c1F_d = c1F_d * c1;

%Do we accept the QP solution ?
fx = f(x);
cx = c(x);
Fx = merit(fx, cx, rho);
F0 = merit(fx0, cx0, rho);
if Fx < F0
    nbcall = nbcall + 1;
    dx = d;
    return%the QP solution is accepted
end


s = 0.5;
dx = s*d;
x = x0 + dx;
fx = f(x);
cx = c(x);
Fx = merit(fx, cx, rho);

k = 0;%rank of the iteration
expr = Fx - F0 - s*c1F_d;
while (k < nbiter) && (expr >= 0)
    k = k + 1;
    s = s/2;
    dx = s*d;
    x = x0 + dx;
    fx = f(x);
    cx = c(x);
    Fx = merit(fx, cx, rho);
    expr = Fx - F0 - s*c1F_d;
end

nbcall = nbcall + 2 + k;

if expr >= 0
    warning("The Armijo's condition is not satisfied")
end

end