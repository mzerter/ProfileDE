function fnval = chemo0_fun(times, y, logp)
%  CHEMO_FUN to be called in Chemostat demo file.  There parameter
%  values are supplied in log-scale, but state values in original scale.
p = exp(logp);
fnval = y;
Q = p(3).*y(:,2) + p(4).*y(:,3);
fnval(1) = p(6).*(p(5) - y(:,1)) - p(12).*y(:,2).*y(:,1)./ ...
    (p(10) + y(:,1)) - p(12).*y(:,3).*y(:,1)./(p(11) + y(:,1));
fnval(2) = y(:,2).*(p(9).*p(12).*y(:,1)./ ...
    (p(10) + y(:,1)) - p(3).*p(13).*...
    (y(:,4) + y(:,5))./(p(15) + Q) - p(6));
fnval(3) = y(:,3).*(p(9).*p(12).*y(:,1)./ ...
    (p(11) + y(:,1)) - p(4).*p(13).*(y(:,4) + y(:,5))./ ...
    (p(15) + Q) - p(6));
fnval(4) = y(:,4).*(p(14).*p(13).*Q./ ...
    (p(15) + Q) - (p(6) + p(7) + p(8)));
fnval(5) = p(8).*y(:,4) - (p(6) + p(7)).*y(:,5);

end

