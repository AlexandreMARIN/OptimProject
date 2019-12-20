function [x, lambda, k, normGradLag, fx, cx] = SQP(x0, problem, epsilon, x_inf, x_sup, maxnbiter, maxnbcall, mindx, mindf, getMsgWarn)
%the SQP algorithm
%Input:
%x0 : first estimate of a minimiser, a column vector
%problem : a handle to a function which returns the value of the objective
%and the values of the constraint
%epsilon : tolerance
%x_inf : lower bounds for each coordinate of x (see below what is x)
%x_sup : upper bounds for each coordinate of x
%maxnbiter : maximum number of iterations
%maxnbcall : maximum number of calls for f and c
%mindx : minimum move for x, expressed in Euclidian norm
%mindf : minimum decrease for f(x)
%getMsgWarn : an optional vector with two Boolean values:
%  - if the first value is true, any message can be written
%  - if the second value is true, any warning can be raised
%By default, that last argument is set to [false, false]
%Output:
%x : an estimate of a minimiser
%lambda : multipliers of Lagrange
%k : number of iterations
%normGradLag : Euclidian norm of the gradient of the Lagrangian function
% at x
%fx, cx : the value of the objective and of the constraints at x

if nargin < 10
    getMsgWarn = [false, false];
end


x = x0;

[n, ~] = size(x);
h = 1e-4 * x;%for computing gradients
c1 = 0.1;%see Armijo
nbit_Armijo = 10;%number of iterations for Armijo

k = 0;%rank of the current iteration
nbcall = 1;%number of calls for f, c

%compute A, Q, g, b
H = eye(n);
[fx, cx] = problem(x);
[g, A] = gradient(x, fx, cx, h, problem);
b = -cx;
nbcall = nbcall + n;

%we solve the quadratic problem
[d_QP, lambda] = KKT_quad_lin(H, g, A, b, getMsgWarn(2));
new_x = x + d_QP;

%we make some projections
ind_inf = new_x < x_inf;
new_x(ind_inf) = x_inf(ind_inf);
d_QP(ind_inf) = x_inf(ind_inf) - x(ind_inf);
ind_sup = new_x > x_sup;
new_x(ind_sup) = x_sup(ind_sup);
d_QP(ind_sup) = x_sup(ind_sup) - x(ind_sup);

%globalisation
fprev = fx;%value of f at prev_x
[x, fx, cx, nbcall, dx] = Armijo(d_QP, new_x, x, fx, g, cx, c1, problem, norm(lambda, Inf)+1, nbit_Armijo, nbcall, getMsgWarn(2));


while (k < maxnbiter)

    %compute A, Q, g, b
    h = 1e-4 * x;
    prev_g = g;%previous value of the gradient of the objective
    prev_A = A;%previous value of A
    [g, A] = gradient(x, fx, cx, h, problem);
    gxlag = g + A'*lambda;%gradient of the Lagrangian with respect to x
    nbcall = nbcall + n;
    
    %main stop criterion
    normGradLag = norm([gxlag ; cx]);
    if normGradLag < epsilon
        return
    end
    
    %another three stop criteria
    if nbcall > maxnbcall
        if getMsgWarn(1)
            fprintf("The maximum number of calls for the objective has been reached : \n\tSQP stops\n")
        end
        return
    end
    
    ndx = norm(dx);
    if (ndx > 0) && (ndx < mindx)
        if getMsgWarn(1)
            fprintf("The minimum move for the minimiser has been reached :\n\tSQP stops\n")
        end
        return
    end

    df = fprev - fx;
    if (df > 0) && (df <= mindf)
        if getMsgWarn(1)
            fprintf("The value of the objective function has not enough decreased :\n\tSQP stops\n")
        end
        return
    end
    
    k = k + 1;
    
    b = -cx;
    H = compute_H(H, dx, gxlag, prev_g, prev_A, lambda);

    %we solve the quadratic problem
    [d_QP, lambda] = KKT_quad_lin(H, g, A, b, getMsgWarn(2));
    new_x = x + d_QP;

    %we make some projections
    ind_inf = new_x < x_inf;
    new_x(ind_inf) = x_inf(ind_inf);
    d_QP(ind_inf) = x_inf(ind_inf) - x(ind_inf);
    ind_sup = new_x > x_sup;
    new_x(ind_sup) = x_sup(ind_sup);
    d_QP(ind_sup) = x_sup(ind_sup) - x(ind_sup);
    
    %globalisation
    fprev = fx;
    [x, fx, cx, nbcall, dx] = Armijo(d_QP, new_x, x, fx, g, cx, c1, problem, norm(lambda, Inf)+1, nbit_Armijo, nbcall, getMsgWarn(2));

end

end