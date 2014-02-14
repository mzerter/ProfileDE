function dfdxval = SEIR_dfdx(t,y,p,more)

% These lines for log parameter estimation
% p(1:5) = exp(p(1:5));
beta   = more.beta_fun(t,p,more);
beta   = beta(:);
nvar   = size(y,2);
betaYI = beta.*(p(2) + y(:,3));
betaYS = beta.*y(:,1); 
r = zeros(length(t),nvar,nvar);
r(:,1,1) = -betaYI - p(3);
r(:,1,3) = -betaYS;
r(:,2,1) =  betaYI;
r(:,2,2) = -p(4) - p(3);
r(:,2,3) =  betaYS;
r(:,3,2) =  p(4);
r(:,3,3) = -p(5) - p(3);
dfdxval = r;

end

