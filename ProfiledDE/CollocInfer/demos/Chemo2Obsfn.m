function [x] = Chemo2Obsfn(t,y,p,more)
    y = exp(y);
    
    x = [p(1)*(y(:,2)+y(:,3)), p(2)*(y(:,4)+y(:,5))];
        
    x = log(x);


end