function x = piecewise(t,t1,t2,t3,v1,v2,v3,v4 )
% piecewise : Implement a proper function here!

if (t<=t1)
    x = v1;
else if (t<=t2)
        x = v2;
    else if (t<=t3)
            x = v3;
        else x = v4;
        end;
    end;
end;

end

