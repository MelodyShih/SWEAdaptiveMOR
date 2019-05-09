close all;
clear;
global limiter Np K B

%% Setup solver
% Test the gradient
Globals1D;

% Final time
FinalTime = 10;

% given parameter 
a=0; b=500; % space 

% Polynomial order used for approximation 
N = 1;

scale = 2;
% number of elements 
K = scale*b;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(a, b, K);

limiter = 1; % if use limiter


% Initialize solver and construct grid and metric
StartUp1D;

%% Setup variables for reduced order based
w = 5; % window of size
winit = 2;
wtotal = 2;
n = 10; % number of reduced basis
z = 5;  % how often we adapted the sample pts, set to 1 for testing how well reduced space approximates true solution
m = 2*Np*K; % number of sample points
r = 1;  % rank r updates of the reduced basis
Q = zeros(2*Np*K, wtotal);
F = zeros(2*Np*K, wtotal);

%% Initial condition
% height
hinit = ones(size(x)); % size(x) = (#order of poly+1)x(#discretization pts)

% setup bathymetry (mu)
a0 = (a+b)/2-10;
b0 = (a+b)/2+10;
mu = 0.15;
p = b0 - a0;
B = -hinit+mu*(1 + cos(2*pi/p*(x - (a0+b0)/2))).*(x>a0 & x<b0);

% momentum 
vinit = zeros(Np,K);

%time 
time=0;
% fix time step
CFL=0.1; g=9.8;
mindx = min(abs(x(2,:)-x(1,:)));
tstep = CFL*min(min(mindx./(abs(vinit./hinit)+sqrt(g*hinit))));

%% Plot initial condition
figure(1);
plot(x,hinit+B,'b',x,B,'k','LineWidth',2);
title(['t=',num2str(time)]);

%% Precompute full solution for projection onto reduced space
[Qhfull,Qvfull,time] = solveFOM(hinit,vinit,time,tstep,wtotal);
Qfull = zeros(2*Np*K, wtotal);
Qfull(:,1:wtotal) = [Qhfull;Qvfull];

q = Qfull(:,1);
h = reshape(q(1:Np*K), [Np, K]);
v = reshape(q(Np*K+1:end), [Np, K]);

%%
idx = [1000; 1050];
[Krd,crd,crdm] = ReducePoint(idx');
[F] = StateFixTSRD1D(h,v,tstep,2*tstep,B,limiter,Krd,crd,crdm,idx');

figure;
plot(Qfull(idx,2), 'o');
hold on;
plot(F, '*');
legend('true', 'new');

norm(Qfull(idx,2)-F)

%%
q = Qfull(:,1);
h = reshape(q(1:Np*K), [Np, K]);
v = reshape(q(Np*K+1:end), [Np, K]);
[rhsh, rhsv] = StateRHS1D(h, v, time, B);
rhsh1=rhsh(:,crdm);
rhsv1=rhsv(:,crdm);
q = Qfull(:,1);
h = reshape(q(1:Np*K), [Np, K]);
v = reshape(q(Np*K+1:end), [Np, K]);
[rhsh2, rhsv2] = StateRHSRD1D(h, v, time, B,Krd,crd,crdm,idx');


