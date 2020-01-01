function [x, fx, cx, nbcall, dx, rho] = Armijo(d, x, x0, fx0, gx0, cx0, c1, problem, rho0, x_inf, x_sup, nbiter, nbcall, maxnbcall, getWarn)
%Input:
%d : d = x - x0, direction of descent
%x : the approximation of the minimiser given by the Newton's algorithm
%x0 : the former approximation of the minimiser
%fx0, gx0, cx0 : resp. values of f, gradient of f and value of c at x0
%c1 : for the Armijo's algorithm
%problem : it must be the same thing as in the SQP function
%rho0 : for the merit function, initial value of rho
%x_inf, x_sup : lower and upper bounds for each coordinate of the minimiser
%nbiter : maximum number of iterations in Armijo()
%nbcall : current number of calls for f, c
%maxnbcall : maximum number of calls for f, c
%getWarn : a Boolean value, true means any warning can be send
%
%Output:
%x : the new estimate of the minimiser
%fx, cx : f(x), c(x)
%nbcall : number of calls for f, c
%dx : = x - x0
%rho : value used for the merit function

%we compute the derivative of the merit function,
%with respect to the direction d
%we multiply it by c1
rho = rho0;
normcx0 = norm(cx0, 1)+1;
gx0d = gx0'*d;
c1F_d = gx0d - rho*normcx0;
%we must have a direction of descent
%so we increase the value of rho
while (c1F_d >= 0)
    rho = rho * 10;
    c1F_d = gx0d - rho*normcx0;
end
c1F_d = c1F_d * c1;

%we make some projections
ind_inf = x < x_inf;
x(ind_inf) = x_inf(ind_inf);
ind_sup = x > x_sup;
x(ind_sup) = x_sup(ind_sup);
dx = d;
dx(ind_inf) = x_inf(ind_inf) - x0(ind_inf);
dx(ind_sup) = x_sup(ind_sup) - x0(ind_sup);

%Do we accept the QP solution ?
[fx, cx] = problem{1}(x);
cx = cx(problem{2});
fx_newton = fx;
cx_newton = cx;
x_newton = x;
d_newton = dx;
Fx = merit(fx, cx, rho);
F0 = merit(fx0, cx0, rho);
nbcall = nbcall + 1;
if (Fx < F0) || (nbcall >= maxnbcall)
    return%the QP solution is accepted
end


s = 0.5;
dx = s*d;
x = x0 + dx;

%we make some projections
ind_inf = x < x_inf;
x(ind_inf) = x_inf(ind_inf);
dx(ind_inf) = x_inf(ind_inf) - x0(ind_inf);
ind_sup = x > x_sup;
x(ind_sup) = x_sup(ind_sup);
dx(ind_sup) = x_sup(ind_sup) - x0(ind_sup);

[fx, cx] = problem{1}(x);
cx = cx(problem{2});
Fx = merit(fx, cx, rho);
nbcall = nbcall + 1;

k = 0;%rank of the iteration
expr = Fx - F0 - s*c1F_d;
while (k < nbiter) && (nbcall < maxnbcall) && (expr >= 0)
    k = k + 1;
    s = s/2;
    dx = s*d;
    x = x0 + dx;
    %we make some projections
    ind_inf = x < x_inf;
    x(ind_inf) = x_inf(ind_inf);
    dx(ind_inf) = x_inf(ind_inf) - x0(ind_inf);
    ind_sup = x > x_sup;
    x(ind_sup) = x_sup(ind_sup);
    dx(ind_sup) = x_sup(ind_sup) - x0(ind_sup);

    [fx, cx] = problem{1}(x);
    cx = cx(problem{2});
    Fx = merit(fx, cx, rho);
    expr = Fx - F0 - s*c1F_d;
    nbcall = nbcall + 1;
end


if expr >= 0
    if getWarn
        warning("The Armijo's condition is not satisfied")
    end
    %we keep the estimate given by the Newton's method
    %(after making projections)
    x = x_newton;
    fx = fx_newton;
    cx = cx_newton;
    dx = d_newton;
end

end