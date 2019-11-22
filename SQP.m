function [x, lambda, k, normGradLag] = SQP(x0, f, c, epsilon, x_inf, x_sup, maxnbiter, maxnbcall, mindx, mindf)
%the SQP algorithm
%Input:
%x0 : first estimate of a minimiser
%f, c : objective, constraints
%epsilon : tolerance
%x_inf : lower bounds for each coordinate of x (see below what is x)
%x_sup : upper bounds for each coordinate of x
%maxnbiter : maximum number of iterations
%maxnbcall : maximum number of calls for f and c
%mindx : minimum move for x, expressed in Euclidian norm
%mindf : minimum decrease for f(x)
%Output:
%x : an estimate of a minimiser
%lambda : multipliers of Lagrange
%k : number of iterations
%normGradLag : Euclidian norm of the gradient of the Lagrangian function
% at x

x = x0;

[n, ~] = size(x);
h = 1e-3;%for computing gradients
c1 = 0.1;%see Armijo
nbit_Armijo = 1000;%number of iterations for Armijo

k = 0;%rank of the current iteration
nbcall = 1;%number of calls for f, c

%compute A, Q, g, b
H = eye(n);
fx = f(x);
cx = c(x);
[g, A] = gradient(x, fx, cx, h, f, c);
b = -cx;
nbcall = nbcall + n;

%we solve the quadratic problem
[d_QP, lambda] = KKT_quad_lin(H, g, A, b);
new_x = x + d_QP;

%we make some projections here
ind_inf = new_x < x_inf;
new_x(ind_inf) = x_inf(ind_inf);
ind_sup = new_x > x_sup;
new_x(ind_sup) = x_sup(ind_sup);

%globalisation
prev_x = x;%previous value of x
fprev = fx;%value of f at prev_x
cprev = cx;%value of c at prev_x
[x, fx, cx, nbcall, dx] = Armijo(d_QP, new_x, x, fx, g, cx, c1, f, c, norm(lambda, Inf), nbit_Armijo, nbcall);


while (k < maxnbiter)

    %compute A, Q, g, b
    [g, A] = gradient(x, fx, cx, h, f, c);
    nbcall = nbcall + n;
    
    %main stop criterion
    normGradLag = norm([g + A'*lambda ; cx]);
    if normGradLag < epsilon
        return
    end
    
    %another three stop criteria
    if nbcall > maxnbcall
        fprintf("The maximum number of calls for the objective has been reached : \n\tSQP stops\n")
        return
    end
    
    if norm(dx) < mindx
        fprintf("The minimum move for the minimiser has been reached :\n\tSQP stops\n")
        return
    end

    df = fprev - fx;
    if (df >= 0) && (df <= mindf)
        fprintf("The value of the objective function has not enough decreased :\n\tSQP stops\n")
        return
    end
    
    k = k + 1;
    
    b = -cx;
    H = compute_H(H, prev_x, fprev, cprev, x, fx, cx, lambda, f, c);

    %we solve the quadratic problem
    [d_QP, lambda] = KKT_quad_lin(H, g, A, b);
    new_x = x + d_QP;

    %we make some projections
    ind_inf = new_x < x_inf;
    new_x(ind_inf) = x_inf(ind_inf);
    ind_sup = new_x > x_sup;
    new_x(ind_sup) = x_sup(ind_sup);
    
    %globalisation
    prev_x = x;
    fprev = fx;
    cprev = cx;
    [x, fx, cx, nbcall, dx] = Armijo(d_QP, new_x, x, fx, g, cx, c1, f, c, norm(lambda, Inf), nbit_Armijo, nbcall);

end

end