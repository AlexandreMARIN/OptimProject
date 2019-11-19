function [x, fx, cx, Fx] = Armijo(x, x0, F0, gx0, cx0, c1, f, c, rho, nbiter)
%Input:
%x : the approximation of the minimiser given by the Newton's algorithm
%x0 : the former approximation of the minimiser
%F0 : the evaluation of the merit function at x0
%gx0, cx0 : resp. values of c and gradient of f at x0
%c1 : for the Armijo's algorithm
%f, c : objective function and constraints
%rho : for the merit function
%nbiter : number of iterations
%
%Output:
%x : the new estimate of the minimiser
%fx, cx : f(x), c(x)
%Fx : the value of the merit function at x

d = x - x0;

%Do we accept the QP solution ?
fx = f(x);
cx = c(x);
Fx = merit(fx, cx, rho);
if Fx < F0
    return
end

%the derivative of the merit function, with respect to the direction d
%we multiply it by c1
c1F_d = c1*(gx0'*d - rho*norm(cx0, 1));

s = 0.5;
x = x0 + s*d;
fx = f(x);
cx = c(x);
Fx = merit(fx, cx, rho);

k = 0;%rank of the iteration
while k < nbiter && Fx >= F0 + s*c1F_d
    k = k + 1;
    s = s/2;
    x = x0 + s*d;
    fx = f(x);
    cx = c(x);
    Fx = merit(fx, cx, rho);
end

end