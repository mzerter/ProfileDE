function y =  RMobsfun(t,x,p,more)
    x = exp(x);
    y = [x(1,:) + x(2,:) , x(3,:)];

    y = log(y);

end

