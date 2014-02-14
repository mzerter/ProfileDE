function dfdx = NS_dfdx(times, x, p, more)
betabasis = more.betabasis;
betanbasis = getnbasis(betabasis);
betamat = eval_basis(times, betabasis);
betacoef = p(1:betanbasis);
betavec = betamat * betacoef;
nobs = length(times);
dfdx = zeros(nobs, 1, 1);
dfdx(:,1,1) = -betavec;
end
