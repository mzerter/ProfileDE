function [coefs, iterhist] = DElsqnonlin(coefs, lsopts,  ...
          basis_cell,Ycell,Tcell,wts,lambda,fn,alg,pars,fn_extras)

%  Last modified 18 June 2007

%  Compute initial function and gradient values

[f,J] = SplineCoefErr(coefs,basis_cell,Ycell,Tcell,wts,lambda, ...
                                fn,alg,pars,fn_extras);
SSE  = sum(f.^2);
grad = 2.*(J'*f);
gradnorm = full(sqrt(sum(grad.^2)));

%  compute the initial expected Hessian

hessmat = 2.*(J'*J);

%  evaluate the initial line search direction vector

deltac = -symsolve(hessmat, grad);

%  -------  Initialize iterations  -----------

MAXSTEPITER = 4;
MAXSTEP  = 100;
iternum  = 0;
iterlim  = lsopts.MaxIter;
trial    = 1;
reset    = 0;
linemat  = zeros(3,5);
coefsold = coefs;
SSEold   = SSE;
if strcmp(lsopts.Display, 'off'),    dbglev = 0;  end
if strcmp(lsopts.Display, 'iter'),   dbglev = 1;  end
if strcmp(lsopts.Display, 'notify'), dbglev = 2;  end
status = [iternum, SSE, gradnorm];
iterhist = status;
if dbglev >= 1
    fprintf('%3.f %10.4f %10.4f\n', status(1:3));
end


%  ---------------  beginning of optimization loop  -----------

for iter = 1:iterlim
    iternum = iternum + 1;
    %  initialize logical variables controlling line search
    dblwrd = [0,0];  
    ips    = 0;
    %  compute slope at 0 for line search
    linemat(2,1) = sum(deltac.*grad);
    %  normalize search direction vector
    sdg          = sqrt(sum(deltac.^2));
    deltac       = deltac./sdg;
    linemat(2,1) = linemat(2,1)/sdg;
    % initialize line search vectors
    linemat(:,1:4) = [0; linemat(2,1); SSE]*ones(1,4);
    stepiter  = 0;
    if dbglev >= 2
        fprintf('                 %3.f %10.4f %12.6f %12.6f\n', ...
            [stepiter, linemat(:,1)']);
    end
    %  return with error condition if initial slope is nonnegative
    if linemat(2,1) >= 0
        if dbglev >= 2, disp('Initial slope nonnegative.'); end
        break;
    end
    %  return successfully if initial slope is very small
    if linemat(2,1) >= -lsopts.TolFun;
        if dbglev >= 2, disp('Initial slope too small'); end
        break;
    end
    %  first step set to trial
    linemat(1,5)  = trial;
    %  -------  Begin line search iterations  -----------
    coefsnew = coefs;
    for stepiter = 1:MAXSTEPITER
        if linemat(1,5) <= lsopts.TolX
            %  Current step size too small ... terminate
            if dbglev >= 2
                fprintf('Stepsize too small: %15.7f\n', linemat(1,5));
            end
        end
        %  compute new function value and gradient
        coefsnew = coefs + linemat(1,5).*deltac;  %  update coefficients
        [f,J] = SplineCoefErr(coefsnew,basis_cell,Ycell,Tcell,wts,lambda, ...
                              fn,alg,pars,fn_extras);
        SSE  = sum(f.^2);
        grad = 2.*(J'*f);
        gradnorm = full(sqrt(sum(grad.^2)));
        linemat(3,5) = SSE;
        %  compute new directional derivative
        linemat(2,5) = sum(deltac.*grad);
        if dbglev >= 2
            fprintf('                 %3.f %10.4f %12.6f %12.6f\n', ...
                [stepiter, linemat(:,5)']);
        end
        %  compute next line search step, also testing for convergence
        [linemat, ips, ind, dblwrd] = ...
                              DEstepit(linemat, ips, dblwrd, MAXSTEP);
        trial  = linemat(1,5);
        if trial == MAXSTEP, break; end
        %  ind == 0 means convergence
        if ind == 0 || ind == 5, break; end
    end
    %  -------  End of line search iterations  -----------
    coefs = coefsnew;
    %  check that function value has not increased
    if SSE > SSEold
        %  if it has, terminate iterations with a warning
        if dbglev >= 2
            fprintf('Criterion increased:');
            fprintf('%10.4f %10.4f\n',[SSEold, SSE]);
        end
        %  reset parameters and fit
        coefs  = coefsold;
        SSE    = SSEold;
        deltac = -grad;
        if dbglev > 2
            for i = 1:nbasis, fprintf('%10.4f%', coefs(i)); end
            fprintf('\n');
        end
        if reset == 1
            %  This is the second time in a row that this
            %     has happened ...  quit
            if dbglev >= 2
                fprintf('Reset twice, terminating.\n');
            end
            return;
        else
            reset = 1;
        end
    else
        %  function value has not increased,  check for convergence
        if abs(SSEold-SSE) < lsopts.TolFun
            status = [iternum, SSE, gradnorm];
            iterhist = [iterhist; status];
            if dbglev >= 1
                fprintf('%3.f %10.4f %10.4f\n', status(1:3));
            end
            break;
        end
        %  update old parameter vectors and fit structure
        coefsold = coefs;
        SSEold = SSE;
        %  update the expected Hessian
        hessmat = 2.*(J'*J);
        %  update the line search direction vector
        deltac = -symsolve(hessmat, grad);
        reset = 0;
    end
    %  store iteration status
    status   = [iternum, SSE, gradnorm];
    iterhist = [iterhist; status];
    if dbglev >= 1
        fprintf('%3.f %10.4f %10.4f\n', status);
    end
end

