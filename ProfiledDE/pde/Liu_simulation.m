pdegplot('squareg')
[p,e,t] = initmesh('squareg'); 
[p,e,t] = refinemesh('squareg',p,e,t);
tlist   = linspace(0,1,10); %% time is from 0 to 10 by step 1
c  = 1;
a  = 1;
u0 = 0; %% initial conditions 

%  By default the boundary conditions are of Dirichlet type,
%  with boundary values fixed at 0.

%  first analysis:  f is 1

f = 1;
u1 = parabolic(u0,tlist,'squareb1',p,e,t,c,a,f,u0);

%  plot the results for each time step

map = jet;
colormap('jet')
for i = 1:length(tlist)
    pdeplot(p,e,t,'xydata',u1(:,i),'mesh','off','colormap','hot');
    title(i);
    pause
end 
%  second analysis:  f = sin(2*pi*time/40) set in 
%  M-file timestep.m

u1 = parabolic(u0,tlist,'squareb1',p,e,t,c,a,'timestep',u0);

%  plot the results

for i = 1:length(tlist)
    pdeplot(p,e,t,'xydata',u1(:,i),'mesh','off','colormap','hot');
    title(i);
    pause
end
