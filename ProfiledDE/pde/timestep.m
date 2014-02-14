function f=timestep(p,t,u,time)
if isnan(time)
    f = NaN;
else
    f = sin(2*pi*time/40);
end