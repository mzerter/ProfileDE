function dfdpval = SEIR_dfdp(t,y,p,more)

% These lines for log parameter estimation
% p(1:5) = exp(p(1:5));
t       = t(:);
beta    = more.beta_fun(t,p,more);
beta    = beta(:);
dbetadp = more.beta_dfdp(t,p,more);
tmpmat = diag(y(:,1).*(p(2)+y(:,3)))*dbetadp;
betaYS = beta.*y(:,1);
r = zeros(length(t),3,length(p));
r(:,1,more.beta_ind) = -tmpmat;
r(:,2,more.beta_ind) =  tmpmat;
r(:,1,1) = 1;
r(:,1,2) = -betaYS;
r(:,2,2) =  betaYS;
r(:,1,3) = -y(:,1);
r(:,2,3) = -y(:,2);
r(:,3,3) = -y(:,3);
r(:,2,4) = -y(:,2);
r(:,3,4) =  y(:,2);
r(:,3,5) = -y(:,3);
% These lines for log parameter estimation
% r(:,:,1) = r(:,:,1).*p(1);
% r(:,:,2) = r(:,:,2).*p(2);
% r(:,:,3) = r(:,:,3).*p(3);
% r(:,:,4) = r(:,:,4).*p(4);
% r(:,:,5) = r(:,:,5).*p(5);




dfdpval  = r;

end

