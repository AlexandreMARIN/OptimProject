function [fx, cx] = pathproblem(theta)
%gives the function and the constraints for the SQP algorithm
%for the path problem
%Remark : we normalize to avoid numeric troubles
global R_c V_p

[R, V, ~, ~] = rocketpath(theta);

fx = -norm(V(end,:))/V_p;

cx = [norm(R(end,:))/R_c - 1 ; (V(end,:)/V_p) * (R(end,:)'/R_c)];

end

