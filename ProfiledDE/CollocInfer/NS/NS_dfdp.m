function dfdp = NS_dfdp(times, x, p, more)
betabasis  = more.betabasis;
betanbasis = getnbasis(betabasis);
phimat = eval_basis(times, betabasis);
dfdp1  = -phimat.*repmat(x, 1, betanbasis);
alphabasis  = more.alphabasis;
alphanbasis = getnbasis(alphabasis);
phimat = eval_basis(times, alphabasis);
rain   = eval_fd(times, more.rainfd);
dfdp2  = phimat.*repmat(rain, 1, alphanbasis);
nobs = length(times);
npar = length(p);
dfdp = zeros(nobs, 1, npar);
dfdp(:,1,:) = [dfdp1, dfdp2];
end
