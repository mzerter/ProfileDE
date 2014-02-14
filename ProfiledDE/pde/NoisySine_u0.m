function u0 = noisysine_u0(p,t,u,time) 

u0 = eval_bifd(p(1,:),p(2,:),noisysinebifd)';
