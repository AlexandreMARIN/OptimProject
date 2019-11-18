function H2 = BFGS(H1, y, d)
%it builds a symmetric positive definite matrix H2
%given a n-by-n symmetric positive definite matrix H1
%and two vectors y and d of size n
aux = y'*d;

if aux<=0
    H2 = H1;
    return
end

aux2 = H1*d;
H2 = H1 + y*y'/aux - aux2*aux2'/(d'*H1*d);

end