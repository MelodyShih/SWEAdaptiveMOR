close all;
clear;
global limiter Np K B

%% Setup solver
% Test the gradient
Globals1D;

% Final time
FinalTime = 1;

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
winit = 150;
wtotal = 313; %(1,63), (2,126), (5,313), (10,627)
n = 100; % number of reduced basis
w = n+1; % window of size
z = 1;  % how often we adapted the sample pts, set to 1 for testing how well reduced space approximates true solution
% m = 2*Np*K; % number of sample points
m = 2000;
r = 1;  % rank r up dates of the reduced basis
Q = zeros(2*Np*K, wtotal);
F = zeros(2*Np*K, wtotal);

debug = 1;

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
if(debug)
    tic;
    [Qhfull,Qvfull,~] = solveFOM(hinit,vinit,time,tstep,wtotal+1);
    toc
    Qfull = zeros(2*Np*K, wtotal+1);
    Qfull(:,1:wtotal+1) = [Qhfull;Qvfull];
end

%% Start simulation
time = 0;
[Qh, Qv, time] = solveFOM(hinit, vinit, time, tstep, winit);
Q(:, 1:winit) = [Qh;Qv];
% norm(Qfull(:,1:winit)-Q(:,1:winit)) % Just checking that Qfull is the same as Q in the window

[U, S] = svd(Q(:,1:winit), 'econ');

Uk = U(:, 1:n);
Pk = qdeim(Uk);

F(:,1:winit) = Q(:,1:winit);
qold = Uk'*Q(:,winit);

errs = zeros(wtotal-(winit),1);

%% AADEIM
tic;
for k = winit+1:wtotal
    time_end = time+tstep;
    qnew = ftilde(Q(:,k-1),time,time_end,Uk,Pk); % Ck = (Pk'*Uk)\Pk'*Fk
    Q(:,k) = Uk*qnew;
    
    if(debug)
        if(mod(k, 10) == 0)
            plotsol(x, Q(:,k), Qfull(:,k), time);
        end
        
        fprintf("||UUtQfull(k)-Qfull(k)|| = %e, ", norm(Uk*Uk'*Qfull(:,k) - Qfull(:,k)));
        fprintf("||Qapprox(k)-Qfull(k)||/||Qfull(k)||= %e\n", norm(Q(:,k) - Qfull(:,k))/norm(Qfull(:,k)));
        errs(k-(winit)) = norm(Uk*Uk'*Qfull(:,k) - Qfull(:,k)); % how well the next solution can be represented in the new basis
        if(errs(k-(winit)) )
    end
    
    if (mod(k, z)==0 || k==winit+1)
        fprintf('adapt sample pts ....\n')
        F(:,k) = ftrue(Q(:,k),time,time_end); % F(:, k) is the surrogate of the full model at timestep k+1
        Rk = -F(:,k-w+1:k) + Uk*(Uk(Pk,:)\F(Pk,k-w+1:k));
        [Rksk,sk] = sort(sum((Rk.^2),2),'descend');
        
        % Local coherence
        if(debug)
            figure(3);
            semilogy(Rksk, 'k-.','Linewidth',1.5);
            title(['$m=$', num2str(m), ', time step $k = $', num2str(k)], 'Interpreter', 'latex', 'FontSize',20);
            xlabel('index $j_i$', 'Interpreter', 'latex', 'FontSize',20);
            legend({'$((U_k C_k - F_k)^Te_{j_i})^2$'}, 'Interpreter', 'latex', 'FontSize',20);
            ylim([1e-15, 1]);
            grid on;
        end
        
        skhat = sk(m+1:end);
        sk = sk(1:m);
    else
%         F(sk,k) = ftrueSk(Q(:,k),time,time_end,sk');
        temp = ftrue(Q(:,k),time,time_end);
        F(sk,k) = temp(sk);
        F(skhat,k) = Uk(skhat,:)*pinv(Uk(sk,:))*F(sk,k);
    end
    
    time = time_end;
    
    Fk = F(:, k-w+1:k);
    [Uk, Pk, ~] = adeim(Uk, Pk, sk, Fk(Pk,:), Fk(sk,:), r);
end
toc
%% plot qtrue - UU^Tqtrue
figure(4);
semilogy(1:wtotal-(winit),(errs'),'k.-','Linewidth',1.5);
title("$\|q_{true}-U_kU_k^Tq_{true}\|_2$, number of basis $n = 10$", 'Interpreter', 'latex', 'FontSize',20);
xlabel('time step $k$','Interpreter', 'latex', 'FontSize',20);
grid on;