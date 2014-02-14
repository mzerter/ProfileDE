function dfdxval = SEIR_var_dfdx(t,y,p,more)

p = more.p.fun(t,more.pdef);
beta = more.beta.fun(t,p,more.betadef);
r = zeros(length(t),size(y,2),size(y,2),size(y,2));
% dimnames(r) = list(NULL,colnames(y),colnames(y),colnames(y))
r(:,1,1,1) = beta.*(p(2)+y(:,3)) + p(3);
r(:,1,1,3) = beta.*y(:,1);
r(:,2,2,1) = beta.*(p(2)+y(:,3));
r(:,2,2,2) = p(4) + p(3);
r(:,2,2,3) = beta.*y(:,1);
r(:,3,3,2) = p(4);
r(:,3,3,3) = p(5) + p(3);
r(:,1,2,1) = -beta.*(p(2)+y(:,3));
r(:,1,2,3) = -beta.*y(:,1);
r(:,2,1,:) = r(:,1,2,:);
r(:,2,3,2) = p(4);
r(:,3,2,:) = r(:,2,3,:);
dfdxval = r;

end

