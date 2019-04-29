function [lh, lv,lu, ltime, Nstep] = State1D(h, v, FinalTime, B)
% 
% function [lh, lv,lu, ltime, Nstep] = State1D
% Purpose  : Integrate 1D shallow water equations with bathymetry until FinalTime starting with
%            initial conditions [h, v]
% third-stage SSP-RK scheme
% Input: h -- the initial state at t=0, length of water 
%        v -- the initial state at t=0, hu momentum
%        FinalTime -- the stopping time
%        B -- the bathymetry
% Output:
%        lh -- list of h at different time, 3-dim matrix
%        lv -- list of v at different time, 3-dim matrix
%        lu -- list of u at different time, 3-dim matrix, u=v/h
%        ltime -- list of time recorded
%        Nstep -- number of steps (or number of times h,v,u recorded)
%        obj -- the value of the objective function

Globals1D;

% Parameters
g=9.8;
CFL = 0.2; % how to decide CFL?
time = 0;

% Prepare for adaptive time stepping
mindx = min(abs(x(2,:)-x(1,:))); 
Nstep =0;

% % Limit initial solution
%   h =SlopeLimitN(h); v=SlopeLimitN(v);

u = v./h;
lh(:,:,1)=h;
lv(:,:,1)=v;
lu(:,:,1)=u;
ltime(1)=time;

% outer time step loop 


while(time<FinalTime)
    dt = CFL*min(min(mindx./(abs(v./h)+sqrt(g*h)))); % how to choose?

  if(time+dt>FinalTime)
    dt = FinalTime-time;
  end

  % 3rd order SSP Runge-Kutta
  
  % SSP RK Stage 1.
  [rhsh,rhsv]  = StateRHS1D(h, v, time, B);
  h1  = h  + dt*rhsh;
  v1 = v + dt*rhsv;
% 
%   % Limit fields
%   h1  = SlopeLimitN(h1); v1 = SlopeLimitN(v1); 

  % SSP RK Stage 2.
  [rhsh, rhsv]  = StateRHS1D(h1, v1, time+dt, B);
  h2   = (3*h  + h1  + dt*rhsh )/4;
  v2  = (3*v + v1 + dt*rhsv)/4;
 
%   % Limit fields
%   h2  = SlopeLimitN(h2); v2 = SlopeLimitN(v2);

  % SSP RK Stage 3.
  [rhsh, rhsv]  = StateRHS1D(h2, v2,time+dt/2, B);
  h  = (h  + 2*h2  + 2*dt*rhsh )/3;
  v = (v + 2*v2 + 2*dt*rhsv)/3;

%   % Limit solution
%   h =SlopeLimitN(h); v=SlopeLimitN(v);

  u = v./h;

  % Increment time and adapt timestep
  time = time+dt;
  Nstep = Nstep +1;
  
  % store the info
  lh(:,:,Nstep+1)=h;
  lv(:,:,Nstep+1)=v;
  lu(:,:,Nstep+1)=u;
  ltime(Nstep+1)=time;
  
  
  figure(1);
  plot(x,h+B,'b',x,B,'k','LineWidth',2);
  title(['t=', num2str(time)]);
end
return
