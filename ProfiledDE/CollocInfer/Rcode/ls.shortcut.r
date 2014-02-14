Smooth.LS <- function(fn, data, times, pars, coefs=NULL, basisvals=NULL,
                       lambda,fd.obj=NULL,more=NULL,weights=NULL,
	                     quadrature=NULL, in.meth='nlminb', control.in=list(), eps=1e-6,
                       posproc=FALSE, poslik=FALSE, discrete=FALSE, names=NULL, sparse=FALSE,
                       likfn = make.id(), likmore = NULL)
{
      
    dims = dim(data)

    profile.obj = LS.setup(pars, coefs, fn, basisvals, lambda,fd.obj, more,
                            data, weights, times, quadrature, eps=1e-6, posproc,
                            poslik, discrete, names=names, sparse=sparse,
                            likfn=make.id(),likmore=NULL)

    lik   = profile.obj$lik
    proc  = profile.obj$proc
    coefs = profile.obj$coefs
    data  = profile.obj$data
    times = profile.obj$times

    Ires   = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
    ncoefs = Ires$coefs
    Ires   = Ires$res
    ncoefs = as.matrix(ncoefs)
    if(!is.null(proc$more$names)){ colnames(ncoefs) = proc$more$names }

    if(!is.null(fd.obj)){
      if(length(dims)>2){
        ncoefs = array(ncoefs,c(length(ncoefs)/(dims[2]*dims[3]),dims[2],dims[3]))
      } else{
         ncoefs = array(ncoefs,c(length(ncoefs)/dims[2],dims[2]))
      }
      fd.obj = fd(ncoefs,fd.obj$basis)
      return( list(fd=fd.obj,lik=lik,proc=proc,inner.result=Ires) )
    }
    else{
      return( list(coefs=ncoefs,lik=lik,proc=proc,inner.result=Ires,data=data,
                   times=times) )
    }
}

###############################################################################


Profile.LS <- function(fn,data,times,pars,coefs=NULL,basisvals=NULL,lambda,
                        fd.obj=NULL,more=NULL,weights=NULL,quadrature=NULL,
                        in.meth='nlminb',out.meth='nls',
                        control.in=list(),control.out=list(),eps=1e-6,
                        active=NULL,posproc=FALSE,poslik=FALSE,discrete=FALSE,
                        names=NULL,sparse=FALSE,likfn=make.id(),likmore=NULL)
{
#    browser()
    if(is.null(active)){ active = 1:length(pars) }

    profile.obj = LS.setup(pars=pars,coefs=coefs,fn=fn,basisvals,lambda=lambda,fd.obj,more,
                           data,weights,times,quadrature,eps=1e-6,
                           posproc,poslik,discrete,names,sparse,
                           likfn=make.id(),likmore=NULL)

    dims = dim(data)

    lik   = profile.obj$lik
    proc  = profile.obj$proc
    coefs = profile.obj$coefs
    data  = profile.obj$data
    times = profile.obj$times
   
    if(file.exists('curcoefs.tmp')){file.remove('curcoefs.tmp')}
    if(file.exists('optcoefs.tmp')){file.remove('optcoefs.tmp')}
    if(file.exists('counter.tmp')){file.remove('counter.tmp')}
    
    Ires = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
    ncoefs = Ires$coefs

    write.table(ncoefs,file='optcoefs.tmp',col.names=FALSE,row.names=FALSE)
    write.table(ncoefs,file='curcoefs.tmp',col.names=FALSE,row.names=FALSE)
    
    apars = pars[active]
    aparamnames = names(apars)
    
    
    if(out.meth == "ProfileGN"){
    ################  Gauss-Newton optimization  #########################
      res=Profile.GausNewt(pars=pars, times=times, data=data, coefs=ncoefs,
		               lik=lik, proc=proc, in.meth=in.meth,
                           control.in=control.in,
		               active=active, control=control.out)
      apars = res$pars[active]

      ncoefs = res$in.res$coefs
      g = res$in.res$df
      resid = res$in.res$f
    }
    else if(out.meth == "nls"){      
    ################  nls optimization  ########################
      if(is.null(control.out$trace)){control.out$trace=TRUE}
      if(is.null(control.out$maxiter)){control.out$maxiter=100}
      if(is.null(control.out$tol)){control.out$tol=1e-8}
      if(is.null(control.out$printEval)){control.out$printEval=TRUE}
      if(is.null(control.out$warnOnly)){control.out$warnOnly=TRUE}
      res = nls(~ProfileSSE(pars, allpars, times, data, coefs, lik, proc,
                            in.meth, control.in, active),
        data = list(allpars=pars, times=times, data=data, coefs=ncoefs,
                    lik=lik, proc=proc,
        in.meth=in.meth,control.in=control.in,active=active),
        start = list(pars=pars[active]),trace=control.out$trace,control=control.out)   
      apars = res$m$getPars()

      g = res$m$gradient()
      resid = res$m$resid()
      if(file.exists('curcoefs.tmp'))
      	 ncoefs = as.matrix(read.table(file='curcoefs.tmp'))
      else 
         ncoefs = coefs
    }
    else{
      res = outeropt(data,times,pars,ncoefs,lik,proc,in.meth=in.meth,out.meth=out.meth,
                      control.in=control.in,control.out=control.out,active=active)
    }

    names(apars) = aparamnames

    pars[active] = apars

    ncoefs = as.matrix(ncoefs)
    if(!is.null(proc$more$names)){ colnames(ncoefs) = proc$more$names }

   
     if(file.exists('curcoefs.tmp')){file.remove('curcoefs.tmp')}
     if(file.exists('optcoefs.tmp')){file.remove('optcoefs.tmp')}
     if(file.exists('counter.tmp')){file.remove('counter.tmp')}
    
    if(!is.null(fd.obj)){
      ncoefs = array(ncoefs,c(nrow(ncoefs)/dims[2],dims[2],dims[3]))
      fd.obj = fd(ncoefs,fd.obj$basis)
      return( list(pars=pars,fd=fd.obj,lik=lik,proc=proc,outer.result=res) )
    }
    else{
      return( list(pars=pars, coefs=ncoefs, lik=lik, proc=proc, outer.result=res,
                   data=data, times=times) )
    }    
}

###############################################################################

LS.setup = function(pars, coefs=NULL, fn, basisvals=NULL, lambda, fd.obj=NULL,
                     more=NULL, data=NULL, weights=NULL, times=NULL,
                     quadrature=NULL, eps=1e-6, posproc=FALSE, poslik = FALSE, 
                     discrete=FALSE, names=NULL, sparse=FALSE,
                     likfn = make.id(), likmore = NULL)
{
	
#  LS_SETUP sets up the LIK, PROC and COEF objects for a CollocInfer
#  analysis where squared error measures of fit to the data and equation
#  are used.   In this common situation, the effort required of the user
#  is considerably reduced by using this function.
#
#  Arguments:  In the follow description of the arguments, N is the
#    number of observations per measured variable, NREP is the number
#    of replications of the observations (usually 1), NBASIS is the
#    number of basis functions used to represent the approximate solution
#    for each variable, and NVAR is the number of variables (including
#    unmeasured variables.)
#
#  PARS       ...  Parameter vector
#  COEFS      ...  A matrix containing coefficients for the expansion. 
#                  As an input argument is either an NBASIS by NVAR matrix,
#                  in which case NREP = 1, or an NBASIS by NVAR by NREP 
#                  array.
#                  On output it has NVAR columns, and its number of rows is
#                  NBASIS*NREP.
#  FN         ...  A struct object containing handles to functions that
#                  evaluate the value of the right side of the 
#                  differential equation as well as the values of a 
#                  number of its partial derivatives.  These functions
#                  must be provided by the user.  However, difference
#                  approximations of derivatives may also be requested.
#  BASISVALS  ...  Either a basis object to be used to approximate each
#                  solution, or a struct object containing the values
#                  of the basis function at each observation time.
#                  Defaults to [].
#  LAMBDA     ...  A single constant, or a vector of length NVAR containing
#                  smoothing parameter values for each variable controlling
#                  how closely a variable comes to solving the equation.
#  FD.OBJ     ...  An alternative means of contributing a basis system.
#                  When this object is supplied, its basis overides that
#                  supplied in argument BASISVALS.  But one or the other
#                  must be supplied.
#                  Defaults to [];
#  MORE       ...  Any additional information required by the functions
#                  referenced in argument FN.
#  DATA       ...  Array of data.  Defaults to [];
#  WEIGHTS    ...  Weights to be applied in computing error sums of 
#                  squares over the observed variables.
#  TIMES      ...  The times of observation of those variables that are
#                  measured.  These times are assumed to be common to all
#                  measured variables.
#  QUADRATURE ...  Locations of quadrature points for approximating the
#                  integrals defining the penalty terms.  If not 
#                  supplied, quadrature points are positioned at the 
#                  centers of the inter-knot intervals.
#  EPS        ...  A small constant value used to compute difference 
#                  approximations to derivatives if this is requested.
#  POSPROC    ...  If nonzero, the functions used in the penalty terms
#                  are constrained to be positive.
#  POSLIK     ...  If nonzero, the functions approximating the data are
#                  constrained to be positive.
#  DISCRETE   ...  If nonzero, a discrete time model is used instead of
#                  the default continuous time model.
#  NAMES      ...  Names for columns
#  SPARSE     ...  Sparse storage model for basis values
#  LIKFUN     ...  Function handle to a defined function for mapping 
#                  trajectories into observations.  Defaults to make_id.
#  LIKMORE    ...  Additional information for LIKFUN


    colnames = names

    if(!is.null(fd.obj)){            # If an fd object is provided, it overrides
      basisvals = fd.obj$basis       # the basis and function values
      
# If an fd object is provided, it overrides values of arguments
#      BASISVALUES and COEFS, if supplied.

     if(!is.null(fd.obj$coefs)){
        coefs = fd.obj$coefs
      }
      if(!is.null(fd.obj$fdnames) & is.null(colnames)){
        colnames = fd.obj$fdnames[[length(fd.obj$fdnames)]]
      }
    }
        
#  ------------------------------------------------------------------------
#  Define the dimensions of the analysis:
#  Determine from coefficient array COEF: 
#      number of basis functions       NBASIS
#      number of repeated measurements NREP
#      number of variables             NVAR
#  Determine from TIMES the number of observations per variable N
#  ------------------------------------------------------------------------

#  set up matrix or array COEFS along with dimensons NBASIS and NREP

   if(length(dim(coefs))>2){
      if(is.null(colnames)){
        colnames = dimnames(coefs)[[3]]
      }
      nrep = dim(coefs)[2]
      #  Put coefs array into matrix format:
      #     Number of rows = NBASIS*NREP
     coefs = matrix(coefs,dim(coefs)[1]*dim(coefs)[2],dim(coefs)[3])
    } else{
        #  COEFS 2-dimensional.  In this case NREP set to 1
        #  and length of second dimension defined as number of variables
      nrep = 1
      if(is.null(colnames)){
        colnames = colnames(coefs)
      }
    }
    
#  If COEFS is empty, set up coefs as zero matrix with NREP = 1

    if(is.null(coefs)){
      if(!is.null(basisvals)){
        if(is.basis(basisvals)) nbasis = basisvals$nbasis 
        else if(is.matrix(basisvals)) nbasis = nrow(basisvals) 
        else nbasis = nrow(basisvals$bvals.obs)
        
        if(!is.null(colnames)) coefs = matrix(0,nbasis,length(colnames))
        else if(!is.null(data)) coefs = matrix(0,nbasis,ncol(data))
        else if(!is.null(weights)) coefs = matrix(0,nbasis,nrow(weights))
        else stop('Cannot determine the dimension of the state vector -- please provide
            one of coefs, fd.obj or names')
      }
    } 

    #  --------------------------------------------------------------------
    #  Define list object LIK containing handles for functions for 
    #  evaluating the values of the error sum of squares for the data and 
    #  their derivatives.
    #  Names of the members of the struct object LIK are:
    #  'fn'      'dfdx'    'dfdy'    'dfdp'    'd2fdx2'  'd2fdxdy'
    #  'd2fdy2'  'd2fdxdp' 'd2fdydp' 'more'    'bvals'
    #  --------------------------------------------------------------------

   #  Level 1 of lik is struct containing functions for evaluating error 
   #  sum of squares for fit to data

   lik = make.SSElik()
    
    #  Add the member MORE to LIK that defines the transformation of the 
    #  process to fit the data.  POS == 0 means no transformation, otherwise an
    #  exponential transformation is applied to provide a positive fit.

    if(!poslik){                        # Map from states to obs
      if(is.list(likfn)){               # All derivatives available analytically
        lik$more = likfn
        if(is.null(lik$more$more)){
          lik$more$more = likmore
        }
      }
      else{                             # Finite-difference for derivatives
        lik$more = make.findif.ode()
        lik$more$more = list(fn=likfn, more=likmore, eps=eps)
      }
    }
    else{                               # States given on log scale
      if(is.list(likfn)){               # Derivatives available analytically
        lik$more = make.exptrans()
        lik$more$more = likfn
        if(is.null(lik$more$more$more)){
          lik$more$more$more = likmore
        }
      }
      else{                             # Finite-difference for derivatives
        lik$more = make.findif.ode()
        lik$more$more$fn = make.exptrans()$fn
        lik$more$more$more = list(fn=likfn,more=likmore, eps=eps)
      }
    }


        
    #  --------------------------------------------------------------------
    #  Define list object PROC containing functions for evaluating the
    #  penalty term and its derivatives
    #  --------------------------------------------------------------------

   #  Level 1 of proc is struct containing functions for evaluating error 
   #  sum of squares for fit to derivative of trajectory

   proc = make.SSEproc()
    
    #  Define a list object PROCMORE containing handles to code for
    #  evalution of right side of differential equation

    if(is.list(fn)){
    #  Argument FN is a struct object with members being function handles 
    #  for functions for evaluating the right side of the
    #  differential equation and its derivatives
      procmore = fn
      procmore$more = more
      findif=FALSE
    }
    else{
    #  Argument FN is a functional data object for evaluating only 
    #  the value of the right side of the differential equation.
    #  Required derivatives are computed by differencing.
      if(is.function(fn)){
        temp.more = list(fn=fn,more=more,eps=eps)
      }
      else if(inherits(fn,'pomp')){
         temp.more = list(fn=pomp.skeleton, more = list(pomp.obj=fn))
      }
      else{
      stop('fn must be either a list of functions, a function or a pomp object')
      }
      
      # base level 1 of procmore contains discrete approximation functions
    
      procmore = make.findif.ode()
      if(posproc){
        #  penalty to be evaluated on log scale, level 2 contains the
        #  logging functions
        procmore$more = list(fn=make.logtrans()$fn,more=temp.more)
      }
      else{
        #  level 2 of procmore contains equation RHS evaluation function
        procmore$more = temp.more
      }
      #  Add an amount defining differences
      procmore$more$eps = eps
      findif=TRUE
    }
    
    
    #  Add  member MORE to PROC containing:
    #    PROCMORE if variables are not constrained to have positive values
    #    the list object returned by function MAKE_LOGTRANS with
    #    an additional member containing PROCMORE

    if(!posproc | findif) {
    #  level 2 of proc is procmore
      proc$more = procmore
    } else { 
      #  level 2 of proc is struct containing logging functions
      proc$more = make.logtrans()
      #  level 3 of proc is procmore
      proc$more$more = procmore
    }
           
    
    proc$more$names = colnames
    proc$more$parnames = names(pars)


#  ------------------------------------------------------------------------
#  Define the basis for representing approximate solutions to the 
#  differential equation, and 
#     add member BVALS to the LIK object
#           to contain basis values and first derivative values
#     define quadrature points and add as member MORE to PROC object
#     define basis values and first derivative
#           values and add as member BVALS to PROC object
#  ------------------------------------------------------------------------

    if(is.basis(basisvals)){
     #  Basis object supplied as an argument.  
     if(is.null(times)){
        stop(paste('if basisvals is is a basis object,',
                   ' you must specify the observation times'))
      }

      #  basis value array for sampling points.  This is always 2-D,
      #  but if there are replications, it is block diagonal with
      #  NREP blocks of basis function values.
      if(sparse){
        lik$bvals = spam(diag(rep(1,nrep)) %x% 
                   eval.basis(times,basisvals))
      } else{
        lik$bvals = diag(rep(1,nrep)) %x% eval.basis(times,basisvals)
      }      
           
      #  define quadrature points
    
      if(is.null(quadrature) | is.null(quadrature$qpts)){
        knots = c(basisvals$rangeval[1],basisvals$params,basisvals$rangeval[2])
        qpts = c(knots[1],knots[-length(knots)] + diff(knots)/2,knots[length(knots)])
      }
      else{
        qpts = quadrature$qpts
      }

      proc$bvals = list()
  
      #  basis value array for quadrature points
    
    if(!discrete){
        #  continuous time model
       if(sparse){
          proc$bvals$bvals  = spam(diag(rep(1,nrep)) %x% 
                                     eval.basis(qpts,basisvals,0))
          proc$bvals$dbvals = spam(diag(rep(1,nrep)) %x%
                                     eval.basis(qpts,basisvals,1))
        }else{
          proc$bvals$bvals  = diag(rep(1,nrep)) %x% eval.basis(qpts,basisvals,0)
          proc$bvals$dbvals = diag(rep(1,nrep)) %x% eval.basis(qpts,basisvals,1)        
        }
        proc$more$weights = matrix(1/length(qpts),length(qpts)*nrep,ncol(coefs))
        #proc$more$weights = 1/length(qpts)
        proc$more$qpts = rep(qpts,nrep)
      }
      else{
        #  discrete time model
       len = length(times)
       basis = eval.basis(times,basisvals,0)
       if(sparse){
         proc$bvals = list(bvals = spam(diag(rep(1,nrep)) %x% basis[1:(len-1),]),
                           dbvals= spam(diag(rep(1,nrep)) %x% basis[2:len,]))
       } else{
         proc$bvals = list(bvals =diag(rep(1,nrep)) %x%  basis[1:(len-1),],
                           dbvals= diag(rep(1,nrep)) %x% basis[2:len,])       
       }
       proc$more$weights = matrix(1/(len-1),(len-1)*nrep,ncol(coefs))
       #proc$more$weights = 1/(len-1)
       proc$more$qpts = rep(times[1:(len-1)],nrep)
       }
      
      
    }
    else{                                   # quadrature is ignored if basisvals
                                            # is not a basis object
     #  Basis values supplied  as argument
     #  quadrature is ignored
     if(discrete & (is.matrix(basisvals) | is.null(basisvals))){
        if(is.null(basisvals)){ basisvals = Diagonal(nrow(coefs)) }
        if(sparse){
          lik$bvals = spam(diag(rep(1,nrep)) %x%basisvals) 
          proc$bvals = list(bvals  = spam(diag(rep(1,nrep)) %x% basisvals[1:(nrow(basisvals)-1),]),
                            dbvals = spam(diag(rep(1,nrep)) %x% basisvals[2:nrow(basisvals),]))        
        }else{
          lik$bvals = diag(rep(1,nrep)) %x%basisvals  
          proc$bvals = list(bvals  = diag(rep(1,nrep)) %x% basisvals[1:(nrow(basisvals)-1),],
                            dbvals =diag(rep(1,nrep)) %x% basisvals[2:nrow(basisvals),])
        }
        proc$more$weights = matrix(1/(nrow(basisvals)-1),
                                   nrow(basisvals)-1,ncol(coefs))
        #proc$more$weights = 1/(nrow(basisvals)-1)
        proc$more$qpts = rep(times[1:(length(times)-1)],nrep)
      }                                    
      else{     
        if(sparse){                                 
          lik$bvals = spam(diag(rep(1,nrep))%x%basisvals$bvals.obs)
    
          proc$bvals =  list(bvals=spam(diag(rep(1,nrep)) %x% 
                                          basisvals$bvals),
                            dbvals=spam(diag(rep(1,nrep)) %x%
                                          basisvals$dbvals))
        } else{
          lik$bvals = diag(rep(1,nrep))%x%basisvals$bvals.obs
    
          proc$bvals =  list(bvals=diag(rep(1,nrep)) %x% basisvals$bvals,
                            dbvals= diag(rep(1,nrep)) %x%basisvals$dbvals)        
        }
        proc$more$qpts = rep(basisvals$qpts,nrep)
        
        if(!is.null(basisvals$qwts))
          proc$more$weights = matrix(basisvals$qwts, 
                                     length(basisvals$qwts)*nrep,ncol(coefs))
        else
          #proc$more$weights = 1/length(proc$more$qpts)
          proc$more$weights = matrix(1/length(proc$more$qpts),
                                       length(proc$more$qpts)*nrep,ncol(coefs))    
      }
    }
 
    #  check consistency of data dimensions with dimensions of other arrays

    if(!is.null(data)){
      if(length(dim(data))==2){
        if(nrep>1){stop('data dimensions must match coefficient dimensions')}
        if(dim(data)[1] != length(times) ){
        stop('data dimensions, times and coefficient dimensions do not match')}
      }
      if(length(dim(data))==3){
         if(dim(data)[2] != nrep | dim(data)[1]!=length(times)){
        stop('data dimensions, times and coefficient dimensions do not match')}
        data = matrix(data,dim(data)[1]*dim(data)[2],dim(data)[3])
        times = rep(times,nrep)
     }
   }
    
   #  ------------------------------------------------------------------------
   #  Define weights for observations and add them in member MORE
   #    of the LIK object
   #  ------------------------------------------------------------------------

   if( is.null(weights) ){
      if(!is.null(data)) lik$more$weights = matrix(1,nrow(data),ncol(data))
      else  lik$more$weights=NULL
   }
   else{
      if(length(dim(weights)) == 3)  weights = matrix(weights,dim(weights)[1]*dim(weights)[2],dim(weights)[3])

      if(!is.null(data)){
        weights =  checkweights(weights,NULL,data)
      }
      else if(length(dim(weights)) == 2)
        lik$more$weights = matrix(weights,dim(weights)[1]*nrep,dim(weights)[2])
      else lik$more$weights = weights
   }
    
  #  ------------------------------------------------------------------------
  #  Define vector LAMBDA of smoothing parameters and add as a
  #  member of member MORE of the PROC object as weights multiplied by 
  #  values in LAMBDA
  #  ------------------------------------------------------------------------

  if(length(lambda) > 1){ proc$more$weights = proc$more$weights%*%diag(lambda) }
   else{ proc$more$weights = as.numeric(lambda)*proc$more$weights }

   return( list(lik=lik,proc=proc,coefs=coefs,data=data,times=times) )

}
