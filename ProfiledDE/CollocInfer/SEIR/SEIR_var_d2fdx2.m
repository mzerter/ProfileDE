function d2fdx2val = SEIR_var_d2fdx2(t,y,p,more)

beta = more.beta.fun(t,p,more.betadef);
r = zeros(length(t),rep(size(y,2),4));
% dimnames(r) = list(NULL,colnames(y),colnames(y),
%         colnames(y),colnames(y))
r(:,1,1,1,3) = beta;
r(:,1,1,3,1) = beta;
r(:,2,2,1,3) = beta;
r(:,2,2,3,1) = beta;
r(:,1,2,1,3) = -beta;
r(:,1,2,3,1) = -beta;
r(:,2,1,:,:) = r(:,1,2,:,:);
d2fdx2val = r;

end

