function [fx, cx] = fc_test1(x)
%for the test case 1
%x must be a real column vector of size 5

fx = (x(1) - 1)^2 + (x(1) - x(2))^2 + (x(2) - x(3))^3 + (x(3) - x(4))^4 + (x(4) - x(5))^4;

cx = [x(1) + x(2)^2 + x(3)^2 - 3*sqrt(2) - 2 ...
    ; x(2) - x(3)^2 + x(4) - 2*sqrt(2) + 2 ...
    ; x(1)*x(5) - 2 ...
    ];
end