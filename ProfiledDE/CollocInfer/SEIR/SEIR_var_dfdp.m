function dfdpval = SEIR_var_dfdp(t,y,p,more)

beta = more.beta_fun(t,p,more);
beta = beta(:);
dbetadp = more.beta_dfdp(t,p,more);
n = size(y,2);
r = zeros(length(t),n,n,length(p));
r(:,1,1,more.beta_ind) =  diag(y(:,1)*(p(2)+y(:,3))) * dbetadp;
r(:,2,2,more.beta_ind) =  diag(y(:,1)*(p(2)+y(:,3))) * dbetadp;
r(:,1,2,more.beta_ind) = -diag(y(:,1)*(p(2)+y(:,3))) * dbetadp;
r(:,2,1,more.beta_ind) = r(:,1,2,more.beta_ind);
r(:,1,1,1) = 1;
r(:,1,1,2) =  beta.*y(:,1);
r(:,1,2,2) = -beta.*y(:,1);
r(:,2,1,2) = r(:,1,2,2);
r(:,2,2,2) =  beta.*y(:,1);
r(:,1,1,3) = y(:,1);
r(:,2,2,3) = y(:,2);
r(:,3,3,3) = y(:,3);
r(:,2,2,4) = y(:,2);
r(:,2,3,4) = y(:,2);
r(:,3,2,4) = y(:,2);
r(:,3,3,5) = y(:,3);
dfdpval = r;

end

