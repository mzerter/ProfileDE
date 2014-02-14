function d2fdxdpval = SEIR_var_d2fdxdp(t,y,p,more)

beta = more.more.beta_fun(t,p,more);
t    = t(:);
beta = beta(:);
dbetadp = more.more.beta_dfdp(t,p,more);
n = size(y,2);
r = zeros(length(t),n,n,n,length(p));
tmpdiag1 = diag(p(2)+y(:,3));
tmpdiag2 = diag(y(:,1));
r(:,1,1,1,more.more.beta_ind) =  tmpdiag1 * dbetadp;
r(:,1,1,3,more.more.beta_ind) =  tmpdiag2 * dbetadp;
r(:,2,2,1,more.more.beta_ind) =  tmpdiag1 * dbetadp;
r(:,2,2,3,more.more.beta_ind) =  tmpdiag2 * dbetadp;
r(:,1,2,1,more.more.beta_ind) = -tmpdiag1 * dbetadp;
r(:,1,2,3,more.more.beta_ind) = -tmpdiag2 * dbetadp;
r(:,1,1,1,2) =  beta;
r(:,2,2,1,2) =  beta;
r(:,1,2,1,2) = -beta;
r(:,1,1,1,3) = 1;
r(:,2,2,2,3) = 1;
r(:,3,3,3,3) = 1;
r(:,2,2,2,4) = 1;
r(:,2,3,2,4) = 1;
r(:,3,2,3,4) = 1;
r(:,3,3,3,5) = 1;
r(:,2,1,:,:) = r(:,1,2,:,:);
d2fdxdpval = r;

end