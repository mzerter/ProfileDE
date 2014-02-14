function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
set(ax,'DataAspectRatio',[3 2 1]);
set(ax,'PlotBoxAspectRatio',[1 1 1]);
set(ax,'XLim',[-1.5 1.5]);
set(ax,'YLim',[-1 1]);
set(ax,'XTick',[ -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1]);
set(ax,'YTickMode','auto');

% Geometry description:
pdepoly([ 
 -0.396, ...
 -0.604, -0.800, -0.995, -0.795, -0.457, -0.399, -0.404, -0.236, -0.191,...
  0.000,  0.204,  0.371,  0.400,  0.289,  0.208,  0.134,  0.048,  0.004,...
  0.008,  0.106,  0.130,  0.004,  0.506,  0.514,  0.820,  1.000,  1.004,...
  1.000,  0.897,  0.800,  0.755,  0.640,  0.800,  0.800,  0.775,  0.608,...
  0.432,  0.408,  0.542,  0.404,  0.330,  0.297,  0.204,  0.004, -0.061,...
 -0.253, -0.391, -0.551, -0.538],...
[ 0.934,...
  0.820,  0.742,  0.604,  0.493,  0.404,  0.338,  0.228,  0.199,  0.338,...
  0.416,  0.538,  0.461,  0.404,  0.346,  0.212,  0.085, -0.004, -0.183,...
 -0.379, -0.600, -0.795, -0.991, -0.971, -0.991, -0.951, -0.648, -0.404,...
 -0.191, -0.126, -0.195, -0.289, -0.204, -0.004,  0.155,  0.236,  0.359,...
  0.412,  0.469,  0.542,  0.677,  0.628,  0.710,  0.677,  0.673,  0.604,...
  0.600,  0.718,  0.722,  0.804], 'P1');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','P1')

% Boundary conditions:
pdetool('changemode',0)
for (i = 1:50)
  pdesetbd(i, 'neu', 1, '0', '0')
end;

% Mesh generation:
setuprop(pde_fig,'Hgrad',1.3);
setuprop(pde_fig,'refinemethod','regular');
pdetool('initmesh')

% PDE coefficients:
%  the 6th arg gives the time interval and sampling rate
%  the 7th arg gives the state of the system at time 0
pdeseteq(2, '1.0', '0.0', '0', '20', '0:.02:2',        ...
         'cos(2*pi*sqrt(x.^2+y.^2)).*sin(atan2(y,x))', ...
         '0.0', '[0 100]')
setuprop(pde_fig, 'currparam', ['1.0'; '0.0'; '0  '; '20 '])

% Solve parameters:
setuprop(pde_fig,' solveparam',...
         str2mat('0','1000','10','pdeadworst',...
                 '0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
%  Flow version
 setuprop(pde_fig,'plotflags',[1 1 1 1 1 1 6 1 0 0 0 41 1 0 0 1 0 1]);
%  No flow version
 setuprop(pde_fig,'plotflags',[1 1 1 1 1 1 6 1 0 0 0 41 1 0 0 0 0 1]);
%  Animation version
%setuprop(pde_fig,'plotflags',[1 1 1 1 1 1 6 1 0 0 1 41 1 0 0 0 0 1]);
 setuprop(pde_fig,'colstring','');
 setuprop(pde_fig,'arrowstring','');
 setuprop(pde_fig,'deformstring','');
 setuprop(pde_fig,'heightstring','');

% Solve PDE:

pdetool('solve')
