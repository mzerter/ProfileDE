function fnval = SIR_dynamic_ode(t,x,p,more)
%  Right side function for SIR differential equation
%  See page 3 of Wiki entry on compartment models in epidemiology
%  To be called by CollocInfer functions
%  x(1) = S,    x(2) = I,     x(3) = R
%  p(1) = beta, p(2) = nu, p(3) = mu
% Value of N passed by more.N

N = more.N;

r = x;
r(1) = -p(1).*x(1).*x(2)/N + p(3)*(N-x(1)) ;
r(2) =  p(1).*x(1).*x(2)/N - (p(2)+p(3)).*x(2);
r(3) =  p(2).*x(2) - p(3).*x(3);
fnval = r;

end

