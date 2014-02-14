function d2fdxdpval = SEIR_d2fdxdp(t,y,p,more)

% These lines for log parameter estimation
% p(1:5) = exp(p(1:5));
beta = more.beta_fun(t,p,more);
t = t(:);
beta = beta(:);
dbetadp = more.beta_dfdp(t,p,more);
r = zeros(length(t),3,3,5);
tmpdiag1 = diag(p(2)+y(:,3));
tmpdiag2 = diag(y(:,1));
r(:,1,1,more.beta_ind) = -tmpdiag1*dbetadp;
r(:,1,3,more.beta_ind) = -tmpdiag2*dbetadp;
r(:,2,1,more.beta_ind) =  tmpdiag1*dbetadp;
r(:,2,3,more.beta_ind) =  tmpdiag2*dbetadp;
r(:,1,1,2) = -beta;
r(:,2,1,2) =  beta;
r(:,1,1,3) = -1;
r(:,2,2,3) = -1;
r(:,3,3,3) = -1;
r(:,2,2,4) = -1;
r(:,3,2,4) =  1;
r(:,3,3,5) = -1;
% These lines for log parameter estimation
% r(:,:,:,1) = r(:,:,:,1).*p(1);
% r(:,:,:,2) = r(:,:,:,2).*p(2);
% r(:,:,:,3) = r(:,:,:,3).*p(3);
% r(:,:,:,4) = r(:,:,:,4).*p(4);
% r(:,:,:,5) = r(:,:,:,5).*p(5);
d2fdxdpval = r;

end
