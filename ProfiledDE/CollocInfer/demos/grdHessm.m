function [grd, Hess] = grdHessm(f,x,toler,ny) 
%Input
%       - f is the function
%       - x is the differentiation point (vector valued)
%       - toler is the tolerance for the relative error
%Output 
%       - grd is the gradient
%       - Hess is the Hessian

%  last modified 12 April 2011 by Jim

indy = 1:ny;
indp = (ny+1):length(x);
yvec = x(indy)';
pvec = x(indp);

%  set default for toler

if nargin < 3 || isempty(toler), toler = 6e-6;  end

f0  = f(0, yvec, pvec);
p   = length(x);
n   = length(f0);
eps = 2.220446e-16;

%  approximate gradient

grd   = zeros(n,p);
if nargout == 2
    dHess = zeros(n,p);
end
for i = 1:p
    idx   = i;
    h     = 0.1;
    j     = 1;
    hp    = zeros(p,1);
    hp(idx) = h;
    yvec = x(indy)';
    pvec = x(indp);
    f0    = f(0,    yvec, pvec)';
    xnew  = x + hp;
    yvec  = xnew(indy)';
    pvec  = xnew(indp);
    fh    = f(0, yvec, pvec)';
    xnew  = x - hp;
    yvec  = xnew(indy)';
    pvec  = xnew(indp);
    fmh   = f(0, yvec, pvec)';
    if nargout == 2
        Dmat1 = zeros(n,12,12);
        Dmat2 = Dmat1;
        Dmat1(:,1,1) = (fh-fmh)/(2*h);
        Dmat2(:,1,1) = (fh-2*f0+fmh)/(h^2);
        rerr1 = ones(n,1);
        rerr2 = ones(n,1);
        while (any(rerr1 > toler) || any(rerr2 > toler)) && (j < 12)
            h = h/2;
            hp(idx) = h;
            xnew  = x + hp;
            yvec  = xnew(indy)';
            pvec  = xnew(indp);
            fh    = f(0, yvec, pvec)';
            xnew  = x - hp;
            yvec  = xnew(indy)';
            pvec  = xnew(indp);
            fmh   = f(0, yvec, pvec)';
            Dmat1(:,j+1,1) = (fh-fmh)/(2*h);
            Dmat2(:,j+1,1) = (fh-2*f0+fmh)/(h^2);
            for k = 1:j
                Dmat1(:,j+1,k+1) = ...
                    Dmat1(:,j+1,k) + (Dmat1(:,j+1,k) - Dmat1(:,j,k))/((4^k)-1);
                Dmat2(:,j+1,k+1) = ...
                    Dmat2(:,j+1,k) + (Dmat2(:,j+1,k) - Dmat2(:,j,k))/((4^k)-1);
            end
            err1  = abs(Dmat1(:,j+1,j+1) - Dmat1(:,j,j));
            err2  = abs(Dmat2(:,j+1,j+1) - Dmat2(:,j,j));
            rerr1 = 2*err1./(abs(Dmat1(:,j+1,j+1)) + abs(Dmat1(:,j,j)) + eps);
            rerr2 = 2*err2./(abs(Dmat2(:,j+1,j+1)) + abs(Dmat2(:,j,j)) + eps);
            j = j+1;
        end
        grdi   = Dmat1(:,j,j);
        dHessi = Dmat2(:,j,j);
    else
        Dmat1 = zeros(n,12,12);
        Dmat1(:,1,1) = (fh-fmh)/(2*h);
        rerr1 = ones(n,1);
        while any(rerr1 > toler) && j < 12
            h = h/2;
            hp(idx) = h;
            xnew  = x + hp;
            yvec  = xnew(indy)';
            pvec  = xnew(indp);
            fh    = f(0, yvec, pvec)';
            xnew  = x - hp;
            yvec  = xnew(indy)';
            pvec  = xnew(indp);
            fmh   = f(0, yvec, pvec)';
            Dmat1(:,j+1,1) = (fh-fmh)/(2*h);
            for k = 1:j
                Dmat1(:,j+1,k+1) = ...
                    Dmat1(:,j+1,k) + (Dmat1(:,j+1,k) - Dmat1(:,j,k))/((4^k)-1);
            end
            err1  = abs(Dmat1(:,j+1,j+1) - Dmat1(:,j,j));
            rerr1 = 2*err1./(abs(Dmat1(:,j+1,j+1)) + abs(Dmat1(:,j,j)) + eps);
            j = j+1;
        end
        grdi = Dmat1(:,j,j);
        dHessi = [];
    end
    grd(:,i)   = grdi;
    if nargout == 2
        dHess(:,i) = dHessi;
    end
end

%  approximate Hessian if required

if nargout == 2
    idx = [kron((1:p),ones(1,p))',repmat((1:p),1,p)'];
    idx = idx((idx(:,1)-idx(:,2)~=0),:);
    lpi = size(idx,1);
    Hess  = zeros(n,p,p);
    for i = 1:lpi
        disp(['Hess calculation, i = ',num2str(i)])
        idx1 = idx(i,1);
        idx2 = idx(i,2);
        h    = 0.1;
        j    = 1;
        hp1  = zeros(p,1);
        hp2  = hp1;
        hp1(idx1) = h;
        hp2(idx2) = h;
        Dmat = zeros(n,12,12);
        Dmat(:,1,1) = ...
            (f(x+hp1+hp2,yvec, pvec)'+f(x-hp1-hp2,yvec, pvec)'- ...
            f(x-hp1+hp2,yvec, pvec)'-f(x+hp1-hp2,yvec, pvec)')/(4*h^2);
        rerr = ones(n,1);
        while any(rerr > toler) && (j < 12)
            h = h/2;
            hp1(idx1) = h;
            hp2(idx2) = h;
            Dmat(:,j+1,1) = ...
                (f(x+hp1+hp2,yvec, pvec)'+f(x-hp1-hp2,yvec, pvec)'- ...
                f(x-hp1+hp2,yvec, pvec)'-f(x+hp1-hp2,yvec, pvec)')/(4*h^2);
            for k = 1:j
                Dmat(:,j+1,k+1) = Dmat(:,j+1,k)+(Dmat(:,j+1,k)-Dmat(:,j,k))/((4^k)-1);
            end
            err  = abs(Dmat(:,j+1,j+1)-Dmat(:,j,j));
            rerr = 2*err./(abs(Dmat(:,j+1,j+1))+abs(Dmat(:,j,j))+eps);
            j = j+1;
        end
        Hessik = Dmat(:,j,j);
        Hess(:,idx1,idx2) = Hessik;
        Hess(:,idx2,idx1) = Hess(:,idx1,idx2);
    end
    for i = 1:p
        Hess(:,i,i) = dHess(:,i);
    end
else
    Hess = [];
end

if n == 1
    grd  = grd';
    Hess = squeeze(Hess);
end



