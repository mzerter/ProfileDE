function fnval = SEIR_var_fun(t,y,p,more)

beta = more.beta.fun(t,p,more.betadef);
r = zeros(length(t),size(y,2),size(y,2));
% dimnames(r) = list(NULL,colnames(y),colnames(y));
r(:,1,1) = p(1) + beta.*y(:,1).*(p(2)+y(:,3)) + p(3).*y(:,1);
r(:,2,2) = beta.*y(:,1).*(p(2)+y(:,3)) + p(4).*y(:,2) + p(3).*y(:,2);
r(:,3,3) = p(4).*y(:,2) + p(5).*y(:,3) + p(3).*y(:,3);
r(:,1,2) = -beta.*y(:,1).*(p(2)+y(:,3));
r(:,2,1) = r(:,1,2);
r(:,2,3) = p(4).*y(:,2);
r(:,3,2) = r(:,2,3);
fnval = r;

end

