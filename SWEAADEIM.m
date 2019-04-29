%% Setup solver
% Test the gradient
Globals1D;

% Final time
FinalTime = 10;

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

%% Setup variables for reduced order based
w = 5;
winit = 10;
k = 5;
% Qh = zeros(Np*K, winit);
% Qv = zeros(Np*K, winit);

% initial condition
h = ones(size(x)); % size(x) = (#order of poly+1)x(#discretization pts)

% setup bathymetry (mu)
a0 = 220;
b0 = 280;
mu = 0.15;
p = b0 - a0;
B = -h+mu*(1 + cos(2*pi/p*(x - (a0+b0)/2))).*(x>a0 & x<b0);
v = zeros(Np,K);

limiter = 1; % if use limiter

%% Plot initial condition
% fix time step
tstep = 1; % time step
time =0;
figure(1);
plot(x,h+B,'b',x,B,'k','LineWidth',2);
title(['t=',num2str(time)]);


%% Start simulation
[Qh, Qv, time] = solveFOM(h, v, B, time, tstep, winit, limiter);
Q = [Qh;Qv];

[U, ~] = svd(Q, 'econ');
Uk = U(:, 1:k);
pk = qdeim(Uk);

F = Q(time-w+1:time);
q = Uk'*Q(:,time);

%%
while time<FinalTime
    if FinalTime-time<tstep
        time_end = FinalTime;
    else
        time_end = time+tstep;
    end    
    [h, v] = StateFixTS1D(h, v, time, time_end, B, limiter);
%     time = time_end;
%     Qh(:,time) = reshape(h, [Np*K, 1]);
%     Qv(:,time) = reshape(v, [Np*K, 1]);
end