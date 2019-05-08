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
winit = 30;
wtotal = 940;
n = 10; % number of reduced basis
z = 5;  % how often we adapted the sample pts, set to 1 for testing how well reduced space approximates true solution
m = 2*Np*K; % number of sample points
r = 1;  % rank r updates of the reduced basis
Q = zeros(2*Np*K, wtotal);
F = zeros(2*Np*K, wtotal);

%% Initial condition (load data from step 500)
% height
load('data/h0.mat');
% momentum
load('data/v0.mat');
%time 
load('data/time0.mat');

% setup bathymetry (mu)
a0 = (a+b)/2-10;
b0 = (a+b)/2+10;
mu = 0.15;
p = b0 - a0;
B = -ones(size(x))+mu*(1 + cos(2*pi/p*(x - (a0+b0)/2))).*(x>a0 & x<b0);


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


%% Start simulation
[Qh, Qv, time] = solveFOM(hinit, vinit, time, tstep, winit);
Q(:, 1:winit) = [Qh;Qv];
% norm(Qfull(:,1:winit)-Q(:,1:winit)) % Just checking that Qfull is the same as Q in the window

[U, S] = svd(Qfull(:,:), 'econ'); % This should only take the svd for the initial window size, right? 
                                        % Changed it to Qfull for now for
                                        % the tests
s = diag(S);
semilogy(diag(S), '-ko');
save('sall.mat','s');
waitforbuttonpress;


Uk = U(:, 1:n);
Pk = qdeim(Uk);

% F(:,1:w-1) = Q(:,winit-w+2:winit);
F(:,1:winit) = Q(:,1:winit);
qold = Uk'*Q(:,winit);

errs = zeros(wtotal-(winit),1);

%% AADEIM
for k = winit+1:wtotal
    time_end = time+tstep;
    qnew = ftilde(Q(:,k-1),time,time_end,Uk,Pk);
    
    Q(:,k) = Uk*qnew;
    
    if(mod(k, 20) == 0)
        plotsol(x, Q(:,k), Qfull(:,k), time);
    end
    
    fprintf("|| UUtQfull(k) - Qfull(k) || = %e, ", norm(Uk*Uk'*Qfull(:,k) - Qfull(:,k))); 
    fprintf("|| Qapprox(k)  - Qfull(k) || = %e\n", norm(Q(:,k) - Qfull(:,k)));
    errs(k-(winit)) = norm(Uk*Uk'*Qfull(:,k) - Qfull(:,k)); % how well the next solution can be represented in the new basis

    if (mod(k, z)==0 || k==winit+1)
        fprintf('adapt sample pts ....\n')
%         F(:,k) = ftrue(Qfull(:,k-1),time,time_end);
        F(:,k) = ftrue(Q(:,k-1),time,time_end);
        Rk = F(:, k-w+1:k) - Uk*(Uk(Pk,:)\F(Pk,k-w+1:k));
        [~,sk] = sort(sum((Rk.^2),2),'descend');
        skhat = sk(m+1:end);
        sk = sk(1:m);
    else
%         temp = ftrue(Qfull(:,k-1),time,time_end);
        F(sk,k) = 
        temp = ftrue(Q(:,k-1),time,time_end);
        F(sk,k) = temp(sk);
        F(skhat,k) = Uk(skhat,:)*pinv(Uk(sk,:))*F(sk,k);
    end
    
    time = time_end;
    
    sk = 1:m; 
    Fk = Qfull(:, k-w+1:k);
%     Fk = F(:, k-w+1:s);
    [Uk, Pk] = adeim(Uk, Pk, sk, Fk(Pk,:), Fk(sk,:), r);
end

%%
semilogy(1:wtotal-(winit),(errs'),'o-','Linewidth',1.5);
title("|q_{true}-UU^Tq_{true}|");
ylabel('error')