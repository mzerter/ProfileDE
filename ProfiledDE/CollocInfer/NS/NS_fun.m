function f = NS_fun(times,y,p,more)
m2 = 0;
betabasis = more.betabasis;
betanbasis = getnbasis(betabasis);
m1 = m2 + 1;
m2 = m2 + betanbasis;
betacoef = p(m1:m2);
betamat = eval_basis(times, betabasis);
betavec = betamat * betacoef;
alphabasis  = more.alphabasis;
alphanbasis = getnbasis(alphabasis);
m1 = m2 + 1;
m2 = m2 + alphanbasis;
alphacoef = p(m1:m2);
alphamat  = eval_basis(times, alphabasis);
alphavec  = alphamat * alphacoef;
rain      = eval_fd(times, more.rainfd);
f = -betavec.*y + alphavec.*rain;
end

