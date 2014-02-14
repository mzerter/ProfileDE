function d2fdxdp = NS_d2fdxdp(times, x, p, more)
betabasis = more.betabasis;
phimat = eval_basis(times, betabasis);
alphanbasis = getnbasis(more.alphabasis);
nobs = length(times);
npar = length(p);
d2fdxdp = zeros(nobs, 1, 1, npar);
d2fdxdp(:,1,1,:) = [-phimat, zeros(nobs, alphanbasis)];
end
