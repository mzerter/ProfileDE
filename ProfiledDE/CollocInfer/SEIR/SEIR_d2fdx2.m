function d2fdx2val = SEIR_d2fdx2(t,y,p,more)

% These lines for log parameter estimation
% p(1:5) = exp(p(1:5));
beta = more.beta_fun(t,p,more);
beta = beta(:);
r = zeros(length(t),3,3,3);
r(:,1,1,3) = -beta;
r(:,1,3,1) = -beta;
r(:,2,1,3) =  beta;
r(:,2,3,1) =  beta;
d2fdx2val = r;

end

