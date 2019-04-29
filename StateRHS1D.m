function [rhsh, rhsv] = StateRHS1D(h, v, time, B)

% function [rhsh, rhsv] = StateRHS1D(h, v, time, B)
% Purpose  : Evaluate RHS flux in 1D shallow  water STATE equations with
% bathymetry 
% Input: h -- the state, length of water 
%        v -- the state, hu
%        time -- the local time, usually used for test
%        B -- the bathymetry
% Output:
%        rhsh -- RHS of h (update in temporal discretization)
%        rhsv -- RHS of v(or hu) (update in temporal discretization)

Globals1D;

% compute maximum velocity for LF flux 
g=9.8; % gravity constant
lm = abs(v./h)+sqrt(g*h);


% Compute fluxes
hf = v; vf = v.^2./h + g*h.^2/2;

% Compute jumps at internal faces
dh = zeros(Nfp*Nfaces,K);  dh(:)  =  h(vmapM)-  h(vmapP); 
dv = zeros(Nfp*Nfaces,K); dv(:)  = v(vmapM)- v(vmapP);
dhf = zeros(Nfp*Nfaces,K); dhf(:)  = hf(vmapM)- hf(vmapP);
dvf = zeros(Nfp*Nfaces,K); dvf(:) = vf(vmapM)- vf(vmapP);
LFc = zeros(Nfp*Nfaces,K); LFc(:) = max(lm(vmapP),lm(vmapM));

% Compute fluxes at interfaces
dhf(:) = nx(:).*dhf(:)/2.0-LFc(:)/2.0.*dh(:); 
dvf(:) = nx(:).*dvf(:)/2.0-LFc(:)/2.0.*dv(:); 

% Boundary conditions (Choose the proper boundary conditions for specific
% problem)

% hin = h(vmapI); hout = h(vmapO); % Neumann boundary condition for h
% vin = 0.0; vout = 0.0; % Dirichlet boundary condition for v

hin = h(vmapO); hout = h(vmapI); % periodic 
vin = v(vmapO); vout = v(vmapI); % periodic

% % specific bounadry condition for test problem

% Set fluxes at inflow/outflow
hfin = vin; vfin = vin.^2./hin + g*hin.^2/2;
lmI=lm(vmapI)/2; nxI=nx(mapI);
dhf(mapI)=nxI*(hf(vmapI)-hfin)/2.0-lmI*(h(vmapI)-hin);  
dvf(mapI)=nxI*(vf(vmapI)-vfin)/2.0-lmI*(v(vmapI)-vin);

hfout = vout; vfout = vout.^2./hout + g*hout.^2/2;
lmO=lm(vmapO)/2; nxO=nx(mapO);
dhf(mapO)=nxO*(hf(vmapO) - hfout)/2.0-lmO*(h(vmapO)- hout);  
dvf(mapO)=nxO*(vf(vmapO)-vfout)/2.0-lmO*(v(vmapO)-vout);


% compute right hand sides of the PDE's
rhsh  = -rx.*(Dr*hf)  + LIFT*(Fscale.*dhf);
rhsv = -rx.*(Dr*(v.^2./h))-g*h.*(rx.*(Dr*(h+B))) + LIFT*(Fscale.*dvf);
return

