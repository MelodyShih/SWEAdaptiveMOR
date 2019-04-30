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

scale = 1;
% number of elements 
K = scale*b;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(a, b, K);

limiter = 1; % if use limiter


% Initialize solver and construct grid and metric
StartUp1D;

%% Setup variables for reduced order based
w = 5; % window of size
winit = 10;
wtotal = 20; 
n = 5; % number of reduced basis
z = 1; % how often we adapted the sample pts, set to 1 for testing how well reduced space approximates true solution
m = 15; % number of sample points
r = 1; % rank r updates of the reduced basis
Q = zeros(2*Np*K, wtotal);
Fk = zeros(2*Np*K, w);

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
time =0;
% fix time step
CFL=0.2; g=9.8;
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


%% Start simulation
[Qh, Qv, time] = solveFOM(hinit, vinit, time, tstep, winit);
Q(:, 1:winit) = [Qh;Qv];

[U, ~] = svd(Q, 'econ');
Uk = U(:, 1:n);
Pk = qdeim(Uk);

Fk(:,1:w-1) = Q(:,winit-w+2:winit);
qold = Uk'*Q(:,winit);

errs = zeros(wtotal-(winit),1);

%% AADEIM
for k = winit+1:wtotal
    time_end = time+tstep;
    qnew = ftilde(qold,time,time_end,Uk,Pk);
    Qnew = Uk*qnew;
    time = time_end;
    qold = qnew;
    if (mod(k, z) == 0 || k==winit+1)
        Fk(:,w) = ftrue(Qnew,time,time_end+tstep);
        Rk = Fk - Uk*(Uk(Pk,:)\Fk(Pk,:));
        [~,sk] = sort(sum((Rk.^2),2),'descend');
        skhat = sk(m+1:end);
        sk = sk(1:m);
    else
        temp = ftrue(Qnew,time,time_end+tstep);
        Fk(sk,k) = temp(sk);
        Fk(skhat,k) = Uk(skhat,:)*pinv(Uk(sk,:))*Fk(sk,k);
    end
    [Uk, Pk] = adeim(Uk, Pk, sk, Fk(Pk,:), Fk(sk,:), r);
    UUtqnew = (Uk*Uk')*Qfull(:,k);
    errs(k-(winit)) = norm(Qfull(:,k)-UUtqnew);
end

plot(1:wtotal-(winit),(errs'));
