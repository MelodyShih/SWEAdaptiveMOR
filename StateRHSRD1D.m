




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



% Compute fluxes


vmapMrd = vmapM(crd);
vmapPrd = vmapP(crd);
% Compute jumps at internal faces
dh = zeros(Nfp*Nfaces,Krd); dh(:)  =  h(vmapMrd)-  h(vmapPrd); 
dv = zeros(Nfp*Nfaces,Krd); dv(:)  = v(vmapMrd)- v(vmapPrd);
dhf = zeros(Nfp*Nfaces,Krd); dhf(:)=dv(:);
dvf = zeros(Nfp*Nfaces,Krd); dvf(:) = v(vmapMrd).^2./h(vmapMrd)+g*h(vmapMrd).^2/2- (v(vmapPrd).^2./h(vmapPrd)+g*h(vmapPrd).^2/2);
LFc = zeros(Nfp*Nfaces,Krd); LFc(:) = max(abs(v(vmapMrd)./h(vmapMrd))+sqrt(g*h(vmapMrd)),abs(v(vmapPrd)./h(vmapPrd))+sqrt(g*h(vmapPrd)));

% Compute fluxes at interfaces
nrd = nx(:);nrd=nrd(crd);
dhf(:) = nrd.*dhf(:)/2.0-LFc(:)/2.0.*dh(:); 
dvf(:) = nrd.*dvf(:)/2.0-LFc(:)/2.0.*dv(:); 

% Boundary conditions (Choose the proper boundary conditions for specific
% problem)

% hin = h(vmapI); hout = h(vmapO); % Neumann boundary condition for h
% vin = 0.0; vout = 0.0; % Dirichlet boundary condition for v

if crdm(1)==1
hin = h(vmapO); vin = v(vmapO); 
else
    hin=h(crd(1)-1);
    vin=v(crd(1)-1);
end

if crdm(end)==K
    hout = h(vmapI); % periodic 
    vout = v(vmapI); % periodic
else
    hout=h(crd(end)+1);
    vout=v(crd(end)+1);
end

% % specific bounadry condition for test problem

%Set fluxes at inflow/outflow

% hfin = vin; vfin = vin.^2./hin + g*hin.^2/2;
% lmI=abs(v(crd(1))./h(crd(1)))+sqrt(g*h(crd(1)))/2; nxI=nx(crd(1));
% dhf(1)=nxI*(v(crd(1))-hfin)/2.0-lmI*(h(crd(1))-hin);  
% dvf(1)=nxI*(v(crd(1)).^2./h(crd(1))+g*h(crd(1)).^2/2-vfin)/2.0-lmI*(v(crd(1))-vin);
% 
% hfout = vout; vfout = vout.^2./hout + g*hout.^2/2;
% lmO=abs(v(crd(end))./h(crd(end)))+sqrt(g*h(crd(end)))/2; nxO=nx(crd(end));
% dhf(end)=nxO*(v(crd(end))-hfout)/2.0-lmO*(h(end)-hout);  
% dvf(end)=nxO*(v(crd(end)).^2./h(crd(end))+g*h(crd(end)).^2/2-vfout)/2.0-lmO*(v(end)-vout);


% compute right hand sides of the PDE's
rhsh  = -rx(:,crdm).*(Dr*v(:,crdm))  + LIFT*(Fscale(:,crdm).*dhf);
rhsv = -rx(:,crdm).*(Dr*(v(:,crdm).^2./h(:,crdm)))-g*h(:,crdm).*(rx(:,crdm).*(Dr*(h(:,crdm)+B(:,crdm)))) + LIFT*(Fscale(:,crdm).*dvf);
return
