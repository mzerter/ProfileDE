classdef ODE < handle
%% ODE Class  
% Idea is to give a symbolic representation of the ODE. 
% Calculates all the necessary functions based on symbolic package.
% TODO:
%           - Add data property to ODE
%           - Graphing data + ode?                
%           - Generate data from ODE, add noise
%           - log data at source function?        
%           - Export functions to file?
%           - Write unit testing


    
    properties(Access=private)
    %Internal variables not exposed to user.
    
        % Display_Symbolic When computeDerivatives() is run, it populates this as a struc. It
        % is used by displayResult() to display the different computed
        % equations in console.
        Display_Symbolic
        % DE_Symbolic Internal copy of symbolic DE set in 
        % constructor.      
        DE_Symbolic
        % paramvector_Symbolic Internal copy of symbolic vector of parameter variables set in 
        % constructor.
        paramvector_Symbolic
        % statevector_Symbolic Internal copy of symbolic vector of state variables set in 
        % constructor.
        statevector_Symbolic
        % time_Symbolic Internal copy of symbolic time variable set in 
        % constructor. 
        time_Symbolic
        % exp_paramvector_bool Boolean value determines if parameter values will be
        % exponentiated in the equations returned by make(). Defaults to 0.
        exp_paramvector_bool
        % f_struc_representation Contains result of computeDerivatives()
        % in a struc for later use by make().
        f_struc_representation

        
    end
    
    methods
        function thisODE = ODE(f,t,x,p)
            syms asim 
            thisODE.DE_Symbolic = asim;
            thisODE.paramvector_Symbolic = asim;
            thisODE.statevector_Symbolic = asim;
            thisODE.time_Symbolic = asim;
            
            thisODE.DE_Symbolic = f;
            thisODE.paramvector_Symbolic = p;
            thisODE.statevector_Symbolic = x;
            thisODE.time_Symbolic = t;
%            Initial Value of exp_paramvector_bool
            thisODE.exp_paramvector_bool = 0;
            
        end
       
        function functionHandle = returnMatlabFunctions(thisODE)
            for i=1:length(thisODE.DE_Symbolic)
                functionHandle{i} = matlabFunction(thisODE.DE_Symbolic(i));
            end
        end
        
        function SetExponentialParameter(thisODE)
            switch thisODE.exp_paramvector_bool
                case {0} 
                    disp(['NOTE: ODE will now exponentiate the parameter vector!'])
                case {1} 
                    disp(['NOTE: ODE will NO LONGER exponentiate the parameter vector!'])
            end
            thisODE.exp_paramvector_bool = not(thisODE.exp_paramvector_bool);
        end
        
        function displayResult(thisODE)
            disp('Value for dfdx')
            disp(thisODE.Display_Symbolic{1})
            disp('Value for dfdp')
            disp(thisODE.Display_Symbolic{2})
            disp('Value for d2fdx2')
            disp(thisODE.Display_Symbolic{3})
            disp('Value for d2fdxdp')
            disp(thisODE.Display_Symbolic{4})
        end
            
        function computeDerivatives(thisODE)
          % Compute functions required for CollocInfer package;

          % Will end up with more than the necessary number of derivatives
          % from these operations.
          X = [thisODE.statevector_Symbolic,thisODE.paramvector_Symbolic];
          
          
          Jf = thisODE.Jacobian2(thisODE.DE_Symbolic,X);
          
          dfdx = Jf(:,1:length(thisODE.statevector_Symbolic));
          thisODE.Display_Symbolic{1} = dfdx;
%         Need some performance testing on this. Kind of convoluted!

          %Need a vector representing the expected parameters!
          syms trick;
          
          fn_out = reshape(thisODE.DE_Symbolic + trick*thisODE.statevector_Symbolic(1),size(thisODE.DE_Symbolic,1)*size(thisODE.DE_Symbolic,2),1);
          fn_hidden(trick,X,thisODE.time_Symbolic) = fn_out;
          
          f.fn = matlabFunction(fn_hidden);
          
          
          dfdx_out = reshape(dfdx + trick*thisODE.statevector_Symbolic(1),size(dfdx,1)*size(dfdx,2),1);          
          dfdx_hidden(trick,X,thisODE.time_Symbolic) = dfdx_out   ;    
          f.dfdx = matlabFunction(dfdx_hidden);
                   
          
%         Something with null length, need to handle error!
          dfdp = Jf(:,end-length(thisODE.paramvector_Symbolic)+1:end);
          thisODE.Display_Symbolic{2} = dfdp;
          
          dfdp_out = reshape(dfdp + trick*thisODE.statevector_Symbolic(1),size(dfdp,1)*size(dfdp,2),1);
          dfdp_hidden(trick,X,thisODE.time_Symbolic) = dfdp_out;
          f.dfdp = matlabFunction(dfdp_hidden);
         
          d2fdx2 = thisODE.Jacobian2(dfdx,thisODE.statevector_Symbolic);
          thisODE.Display_Symbolic{3} = d2fdx2;
          
          d2fdx2_out = reshape(d2fdx2 + trick*thisODE.statevector_Symbolic(1),size(d2fdx2,1)*size(d2fdx2,2),1);
          d2fdx2_hidden(trick,X,thisODE.time_Symbolic) = d2fdx2_out;
          f.d2fdx2 = matlabFunction(d2fdx2_hidden) ;        
          
          d2fdxdp = thisODE.Jacobian2(dfdx,thisODE.paramvector_Symbolic);
          thisODE.Display_Symbolic{4} = d2fdxdp;
          
          d2fdxdp_out = reshape(d2fdxdp + trick*thisODE.statevector_Symbolic(1),size(d2fdxdp,1)*size(d2fdxdp,2),1);
          d2fdxdp_hidden(trick,X,thisODE.time_Symbolic) = d2fdxdp_out;
          f.d2fdxdp = matlabFunction(d2fdxdp_hidden);         
                 
          thisODE.f_struc_representation = f;
          
           disp( ['NOTE: Finished Calculating derivatives for ODE!'] );
            
        end
        
        function plot(thisODE,trange,x_init,p_init,more)
%           Create a figure for plot?
           
            fn = thisODE.make();
           
            fn_new = @(t,y) fn.fn(t,y,p_init,more);
            ode45(fn_new,trange,x_init);
            
        end
        
        function [T,Y] = generateData(thisODE,trange,x_init,p_init,more,sig);
            fn = thisODE.make();
           
            MU = zeros(1,length(x_init));
            SIGMA = sig*eye(length(x_init));
           

            
            fn_new = @(t,y) fn.fn(t,y,p_init,more);
            [T,Y] = ode45(fn_new,trange,x_init);
            
            randNum = mvnrnd(MU,SIGMA,length(T));
            noiseStd = std(randNum);
            signalMean = mean(Y);
            sig2noise = signalMean ./ noiseStd
            Y = Y + randNum;
            
        end
        
        
        
        function f = make(thisODE)
            f.fn = @thisODE.fn;
            f.dfdx = @thisODE.dfdx;
            f.dfdp = @thisODE.dfdp;
            f.d2fdx2 = @thisODE.d2fdx2;
            f.d2fdxdp = @thisODE.d2fdxdp;
            
        end
        

    end
    
    methods (Access=private)
        
       function result = Jacobian2(~,F,y)
           % Replace the symbolic Jacobian. Will compute Jacobian on each
           % of the columns and append each as rows to result. This
           % functions like the Jacobian in the symbolic package pre 2011?.
           
           numCols = size(F,2);
           
           if numCols == 1
               result = jacobian(F,y);
               return; 
           end
           
           J = jacobian(F(:,1),y);
               
           for i = 2 : numCols
               J = cat(1,J,jacobian(F(:,i),y));
           end
           result = J;
           return;
       end
        
       function result = fn(thisODE,t,x,p,more)
            if thisODE.exp_paramvector_bool == 1
                p = exp(p);
            end
%             Need p to be a row vector.
            if isvector(p) == 0
                disp('p is not a vector!')
            end
            p = reshape(p,1,length(p));
            if isvector(x) == 1
                x = reshape(x,1,length(x));
            end

            t = reshape(t,length(t),1);
            xcol = num2cell(x,1);
            pcol = num2cell(p,1);
            
            res = thisODE.f_struc_representation.fn(0,xcol{:},pcol{:},t);
            if isvector(x) == 1
                result = reshape(res,length(x),1);
            else
                result = reshape(res,length(xcol{1}),length(xcol));
            end
           
        end
        
        function result = dfdx(thisODE,t,x,p,more)
            if thisODE.exp_paramvector_bool == 1
                p = exp(p);
            end 
            
            t = reshape(t,length(t),1);
            xcol = num2cell(x,1);
            pcol = num2cell(p,1);
            
            res = thisODE.f_struc_representation.dfdx(0,xcol{:},pcol{:},t);
            result = reshape(res,length(xcol{1}),length(xcol),length(xcol));
           
        end
        
        function result = dfdp(thisODE,t,x,p,more)
            if thisODE.exp_paramvector_bool == 1
                p = exp(p);
            end
            
            t = reshape(t,length(t),1);
            xcol = num2cell(x,1);
            pcol = num2cell(p,1);
            
            res = thisODE.f_struc_representation.dfdp(0,xcol{:},pcol{:},t);
            result = reshape(res,length(xcol{1}),length(xcol),length(pcol));
           
        end
        
        function result = d2fdx2(thisODE,t,x,p,more)
           if thisODE.exp_paramvector_bool == 1
                p = exp(p);
           end
            
           
            t = reshape(t,length(t),1);
            xcol = num2cell(x,1);
            pcol = num2cell(p,1);
            
            res = thisODE.f_struc_representation.d2fdx2(0,xcol{:},pcol{:},t);
            result = reshape(res,length(xcol{1}),length(xcol),length(xcol),length(xcol));
           
        end
        
        function result = d2fdxdp(thisODE,t,x,p,more)
            if thisODE.exp_paramvector_bool == 1
                p = exp(p);
            end
            
            
            t = reshape(t,length(t),1);
            xcol = num2cell(x,1);
            pcol = num2cell(p,1);
            
            res = thisODE.f_struc_representation.d2fdxdp(0,xcol{:},pcol{:},t);
            result = reshape(res,length(xcol{1}),length(xcol),length(xcol),length(pcol));
           
        end
    end
    
end
