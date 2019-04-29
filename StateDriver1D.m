% Test the gradient
Globals1D;

% Final time
FinalTime = 100;

% given parameter 
a=0; b=500; % space 

% Polynomial order used for approximation 
N = 1;

scale = 1;
% number of elements 
K = scale*b;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(a, b, K);


% Initialize solver and construct grid and metric
StartUp1D;


% initial condition
h = ones(size(x));
B = -1+0.5*(x<200).*(x>100).*exp(-(x-150).^2/200);
v = zeros(Np,K);

% % solve the forward equation
% [lh, lv,lu, ltime, Nstep] = State1D(h, v, FinalTime, B);


limiter = 1; % if use limiter

% fix time step 
tstep = 1; % time step
time =0;
figure(1);
plot(x,h+B,'b',x,B,'k','LineWidth',2);
title(['t=',num2str(time)]);
while time<FinalTime
    if FinalTime-time<tstep
        time_end = FinalTime;
    else
        time_end = time+tstep;
    end
    [h, v] = StateFixTS1D(h, v, time, time_end, B,limiter);
    time = time_end;
    figure(1);
    plot(x,h+B,'b',x,B,'k','LineWidth',2);
    title(['t=',num2str(time)]);
end