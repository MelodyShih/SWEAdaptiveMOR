function [h, v] = StateFixTS1D(h, v, T0, FinalTime, B,limiter)
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
% Output:
%        h -- the state at FinalTime, length of water 
%        v -- the state at FinalTime, hu momentum
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
  [rhsh,rhsv]  = StateRHS1D(h, v, time, B);
  h1  = h  + dt*rhsh;
  v1 = v + dt*rhsv;

  if (limiter)
  % Limit fields
  h1  = SlopeLimitN(h1); v1 = SlopeLimitN(v1); 
  end

  % SSP RK Stage 2.
  [rhsh, rhsv]  = StateRHS1D(h1, v1, time+dt, B);
  h2   = (3*h  + h1  + dt*rhsh )/4;
  v2  = (3*v + v1 + dt*rhsv)/4;
 
  if (limiter)
  % Limit fields
  h2  = SlopeLimitN(h2); v2 = SlopeLimitN(v2);
  end
  
  % SSP RK Stage 3.
  [rhsh, rhsv]  = StateRHS1D(h2, v2,time+dt/2, B);
  h  = (h  + 2*h2  + 2*dt*rhsh )/3;
  v = (v + 2*v2 + 2*dt*rhsv)/3;

  if (limiter)
  % Limit solution
  h =SlopeLimitN(h); v=SlopeLimitN(v);
  end


return
