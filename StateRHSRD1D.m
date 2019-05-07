function [rhsh, rhsv] = StateRHSRD1D(h, v, time, B,Krd,crd,crdm,Sk)

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

vmapMrd = vmapM(crd);
vmapPrd = vmapP(crd);
% Compute jumps at internal faces
dh = zeros(Nfp*Nfaces,Krd); dh(:)  =  h(vmapMrd)-  h(vmapPrd); 
dv = zeros(Nfp*Nfaces,Krd); dv(:)  = v(vmapMrd)- v(vmapPrd);
dhf = zeros(Nfp*Nfaces,Krd); dhf(:)  = hf(vmapMrd)- hf(vmapPrd);
dvf = zeros(Nfp*Nfaces,Krd); dvf(:) = vf(vmapMrd)- vf(vmapPrd);
LFc = zeros(Nfp*Nfaces,Krd); LFc(:) = max(lm(vmapPrd),lm(vmapMrd));

% Compute fluxes at interfaces
nrd = nx(:);nrd=nrd(crd);
dhf(:) = nrd.*dhf(:)/2.0-LFc(:)/2.0.*dh(:); 
dvf(:) = nrd.*dvf(:)/2.0-LFc(:)/2.0.*dv(:); 

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
dhf(2*Krd)=nxO*(hf(2*Krd) - hfout)/2.0-lmO*(h(2*Krd)- hout);  
dvf(2*Krd)=nxO*(vf(2*Krd)-vfout)/2.0-lmO*(v(2*Krd)-vout);


% compute right hand sides of the PDE's
rhsh  = -rx(:,crdm).*(Dr*hf(:,crdm))  + LIFT*(Fscale(:,crdm).*dhf);
rhsv = -rx(:,crdm).*(Dr*(v(:,crdm).^2./h(:,crdm)))-g*h(:,crdm).*(rx(:,crdm).*(Dr*(h(:,crdm)+B(:,crdm)))) + LIFT*(Fscale(:,crdm).*dvf);
return
