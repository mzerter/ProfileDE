function basisobj = MakeBasis(range,nbasis,norder,knots,quadvals,dvalue)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function MakeBasis
%
% Makes a B-spline basis with quadrature points.
%
% INPUTS:
%
% range  - the range of the basis
% 
% nbasis - number of basis points
%
% norder - order of the basis functions
%
% knots  - knots for the bspline basis
%
% quadvals -  quadrature points to put as values for the bases
%
% dvalue   -  order of derivative to store quadrature values up to.
%
% OUTPUT:
%
% basisobj    - a B-spline basis with quadrature values attached.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Last modified 18 June 2007

%  check that at least 3 arguments supplied

if nargin < 3
    error('Number of arguments less than 3.');
end

%  Knots set to be equally spaced by default

if nargin < 4
    knots = linspace(range(1), range(2), nbasis-2);
end

%  create basisobj without quadvals or values

basisobj = create_bspline_basis(range,nbasis,norder,knots);

%  add quadrature values if supplied

if nargin >= 5
    basisobj = putquadvals(basisobj,quadvals);
end

%  add values if supplied

if nargin == 6
    values = cell(1,dvalue+1);
    for ivalue=1:(dvalue+1)
        basisvalues    = eval_basis(quadvals(:,1), basisobj, ivalue-1);
        values{ivalue} = basisvalues;
    end
    basisobj = putvalues(basisobj, values);
end