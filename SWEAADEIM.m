clear;
global limiter Np K B

%% Setup solver
% Test the gradient
Globals1D;

% Final time
FinalTime = 10;

% given parameter 
a=0; b=100; % space 

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
wtotal = 11;
n = 5;
z = 5;
m = 3;
r = 1;
Q = zeros(2*Np*K, wtotal);
Fk = zeros(2*Np*K, w);

% initial condition
hinit = ones(size(x)); % size(x) = (#order of poly+1)x(#discretization pts)

% setup bathymetry (mu)
a0 = (a+b)/2-10;
b0 = (a+b)/2+10;
mu = 0.15;
p = b0 - a0;
B = -hinit+mu*(1 + cos(2*pi/p*(x - (a0+b0)/2))).*(x>a0 & x<b0);
vinit = zeros(Np,K);

limiter = 1; % if use limiter

%% Plot initial condition
% fix time step
CFL=0.2; g=9.8;
mindx = min(abs(x(2,:)-x(1,:)));
tstep = CFL*min(min(mindx./(abs(vinit./hinit)+sqrt(g*hinit))));

time =0;
figure(1);
plot(x,hinit+B,'b',x,B,'k','LineWidth',2);
title(['t=',num2str(time)]);


%% Start simulation
[Qh, Qv, time] = solveFOM(hinit, vinit, time, tstep, winit);
Q(:, 1:winit) = [Qh;Qv];

[U, ~] = svd(Q, 'econ');
Uk = U(:, 1:n);
Pk = qdeim(Uk);

Fk(:,1:w-1) = Q(:,winit-w+2:winit);
qold = Uk'*Q(:,winit);

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

end