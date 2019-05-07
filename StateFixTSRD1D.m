function [F] = StateFixTSRD1D(h, v, T0, FinalTime, B,limiter,Krd,crd,crdm,Sk)
% 
% function [lh, lv,lu, ltime, Nstep] = State1D
% Purpose  : Integrate 1D shallow water equations with bathymetry until FinalTime starting with
%            initial conditions [h, v]
% third-stage SSP-RK scheme
% Input: h -- the initial state at t=T0, length of water 
%        v -- the initial state at t=T0, hu momentum
%        T0 -- initial time
%        FinalTime -- final time
%        B -- the bathymetry
%        Sk -- a 2Np*K row vector gives the location of the points
%        Krd -- the number of columns we use 
%        crd -- the column number we use
% Output:
%        F -- the results evluated at Sk
Globals1D;

% Parameters
g=9.8;

time = T0;



if (limiter)
% Limit initial solution
  h =SlopeLimitN(h); v=SlopeLimitN(v);
end

% outer time step loop 

dt=FinalTime-T0;

  % 3rd order SSP Runge-Kutta
  
  % SSP RK Stage 1.
  [rhsh,rhsv]  = StateRHSRD1D(h, v, time, B,Krd,crd,crdm,Sk);
  h1  = h;
  h1(:,crdm)=h1(:,crdm)+ dt*rhsh;
  v1 = v;
  v1(:,crdm)=v1(:,crdm)+ dt*rhsv;

  if (limiter)
  % Limit fields
  h1  = SlopeLimitN(h1); v1 = SlopeLimitN(v1); 
  end

  % SSP RK Stage 2.
  [rhsh, rhsv]  = StateRHSRD1D(h1, v1, time+dt, B,Krd,crd,crdm,Sk);
  h2 = (3*h  + h1)/4;
  h2(:,crdm)=h2(:,crdm)+ dt*rhsh/4;
  v2 = (3*v  + v1)/4;
  v2(:,crdm)=v2(:,crdm)+ dt*rhsv/4;
 
  if (limiter)
  % Limit fields
  h2  = SlopeLimitN(h2); v2 = SlopeLimitN(v2);
  end
  
  % SSP RK Stage 3.
  [rhsh, rhsv]  = StateRHSRD1D(h2, v2,time+dt/2, B,Krd,crd,crdm,Sk);
  h = (h  + 2*h2)/3;
  h(:,crdm) = h(:,crdm)+2*dt*rhsh/3;
  v = (v  + 2*v2)/3;
  v(:,crdm) = v(:,crdm)+2*dt*rhsv/3;


  if (limiter)
  % Limit solution
  h =SlopeLimitN(h); v=SlopeLimitN(v);
  end
  

  F=[h(:);v(:)];
  F = F(Sk);

return
