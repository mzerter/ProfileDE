function fnval = SEIR_ode(t,y,parms)

p    = parms.p;
more = parms.more;

beta = more.beta_fun(t,p,more);
beta = beta(:);
tmpvec = beta.*y(:,1).*(p(2) + y(:,3));
r = y;
r(:,1) = -tmpvec + p(1)              - p(3).*y(:,1);
r(:,2) =  tmpvec - p(4).*y(:,2)      - p(3).*y(:,2);
r(:,3) = p(4).*y(:,2) - p(5).*y(:,3) - p(3).*y(:,3);
fnval = r;

end

