function [x, lambda, k, normGradLag, fx, cx, rho, nbcall] = SQP(x0, problem, epsilon, x_inf, x_sup, maxnbiter, maxnbcall, mindx, mindf, getMsgWarn)
%the SQP algorithm searches for a minimiser x
%Input:
%x0 : first estimate of a minimiser, a column vector
%problem : it is
%   - either a handle to a function 'f' which returns the value of the 
%     objective and the values of the constraint 'c'
%   - or a 2D array cell of the form
%     {handle to a function as above ; indices of the selected constraints}
%  no index in the second cell means "no constraint".
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
%rho : penalty used in Armijo()
%nbcall : number of calls for f and c


%Will we display information ?
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

%we compute A, Q, g, b to build the linear system
H = eye(n);%an estimate of the Hessian matrix Q of the Lagrangian function
%we see what kind of data is 'problem'
%'problem' will be an array cell of the form:
%{handle to the objective-constraints function ; indices of the constraints}
if iscell(problem)
    [fx, cx] = problem{1}(x);
    cx = cx(problem{2});%we take only the constraints given by problem{2}
else
    [fx, cx] = problem(x);
    problem = {problem ; 1:size(cx, 1)};%we take all the constraints
end
%we compute the gradient 'g' of the objective
%and the Jacobian matrix 'A' of the constraints
[g, A] = gradient(x, fx, cx, h, problem);
b = -cx;
nbcall = nbcall + n;%the increase is due to the computation of gradients

%iteration 0 : we compute the multipliers of Lagrange
%with the least squares method
if isempty(cx)
    lambda = [];
    normGradLag = norm(g);%no constraint !
else
    [lambda, ~] = linsolve(A*A', -A*g);
    gxlag = g + A'*lambda;%gradient of the Lagrangian function with respect to x
    %norm of the gradient of the Lagrangian function with respect to (x, lambda)
    normGradLag = norm([gxlag ; cx]);
end

if getMsgWarn(1)
    %we display some information
    fprintf("iteration: %d\n"...
        +"number of calls: %d\n"...
        +"norm of the gradient of the Lagrangian function: %f\n"...
        +"x: ( %s )\n"...
        +"objective at x: %f\n"...
        +"constraints at x: ( %s )\n"...
        +"lambda: ( %s )\n\n"...
        , k...
        , nbcall...
        , normGradLag...
        , join(string(x), ", ")...
        , fx...
        , join(string(cx), ", ")...
        , join(string(lambda), ", "));
end

%the main stop criteria:
if normGradLag < epsilon
    rho = NaN;
    return
end


k = 1;%the first iteration

while (k <= maxnbiter)

    %we solve the quadratic problem
    [d_QP, lambda] = KKT_quad_lin(H, g, A, b, getMsgWarn(2));
    new_x = x + d_QP;

    %globalisation
    fprev = fx;%value of f at prev_x
    [x, fx, cx, nbcall, dx, rho] = Armijo(d_QP, new_x, x, fx, g, cx, c1, problem, norm(lambda, Inf)+1, x_inf, x_sup, nbit_Armijo, nbcall, maxnbcall, getMsgWarn(2));

    %Is the maximum number of calls for 'problem' reached ?
    if nbcall > maxnbcall
        if getMsgWarn(1)
            fprintf("The maximum number of calls for the objective has been reached : \n\tSQP stops\n")
        end
        return
    end

    %compute A, Q, g, b
    h = 1e-4 * x;
    prev_g = g;%previous value of the gradient of the objective
    prev_A = A;%previous value of A
    [g, A] = gradient(x, fx, cx, h, problem);
    
    if isempty(lambda)
        gxlag = g;%no constraint
    else
        gxlag = g + A'*lambda;%gradient of the Lagrangian function with respect to x
    end
    nbcall = nbcall + n;
    
    %for the main stop criterion
    normGradLag = norm([gxlag ; cx]);
    
    if getMsgWarn(1)
        %we display some information
        fprintf("iteration: %d\n"...
               +"number of calls: %d\n"...
               +"norm of the gradient of the Lagrangian function: %f\n"...
               +"x: ( %s )\n"...
               +"objective at x: %f\n"...
               +"constraints at x: ( %s )\n"...
               +"lambda: ( %s )\n"...
               +"penalisation: %f\n\n"...
               , k...
               , nbcall...
               , normGradLag...
               , join(string(x), ", ")...
               , fx...
               , join(string(cx), ", ")...
               , join(string(lambda), ", ")...
               , rho);
    end
    
    if normGradLag < epsilon
        return
    end
    
    %another two stop criteria
    
    %Is the update of x negligible ?
    ndx = norm(dx);
    if (ndx >= 0) && (ndx < mindx)
        if getMsgWarn(1)
            fprintf("The minimum move for the minimiser has been reached :\n\tSQP stops\n")
        end
        return
    end

    %Is the decrease of f(x) negligible ?
    df = fprev - fx;
    if (df >= 0) && (df <= mindf)
        if getMsgWarn(1)
            fprintf("The value of the objective function has not enough decreased :\n\tSQP stops\n")
        end
        return
    end
    
    k = k + 1;

    %we prepare the next iteration
    b = -cx;
    H = compute_H(H, dx, gxlag, prev_g, prev_A, lambda);

end

%necessarily, if we are here the maximum number of iterations is reached
if getMsgWarn(1)
    fprintf("The maximum number of iterations has been reached :\n\tSQP stops\n")
end

end