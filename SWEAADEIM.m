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
w = 5; % window of size
winit = 20;
wtotal = 627; %(1,63), (2,126), (5,313), (10,627)
n = 10; % number of reduced basis
z = 5;  % how often we adapted the sample pts, set to 1 for testing how well reduced space approximates true solution
m = 2*Np*K; % number of sample points
r = 1;  % rank r updates of the reduced basis
Q = zeros(2*Np*K, wtotal);
F = zeros(2*Np*K, wtotal);

debug = 0;

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
tic;
[Qhfull,Qvfull,~] = solveFOM(hinit,vinit,time,tstep,wtotal+1);
toc
Qfull = zeros(2*Np*K, wtotal+1);
Qfull(:,1:wtotal+1) = [Qhfull;Qvfull];

%% Start simulation
time = 0;
[Qh, Qv, time] = solveFOM(hinit, vinit, time, tstep, winit);
Q(:, 1:winit) = [Qh;Qv];
% norm(Qfull(:,1:winit)-Q(:,1:winit)) % Just checking that Qfull is the same as Q in the window

[U, S] = svd(Qfull(:,1:winit), 'econ'); % This should only take the svd for the initial window size, right?
% Changed it to Qfull for now for
% the tests

Uk = U(:, 1:n);
Pk = qdeim(Uk);

% F(:,1:w-1) = Q(:,winit-w+2:winit);
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
    end
    
    if (mod(k, z)==0 || k==winit+1)
        fprintf('adapt sample pts ....\n')
        %         F(:,k) = ftrue(Qfull(:,k-1),time,time_end);
        F(:,k) = ftrue(Q(:,k),time,time_end); % F(:, k) is the surrogate of the full model at timestep k+1
        Rk = -F(:,k-w+1:k) + Uk*(Uk(Pk,:)\F(Pk,k-w+1:k));
        [~,sk] = sort(sum((Rk.^2),2),'descend');
        skhat = sk(m+1:end);
        sk = sk(1:m);
    else
        %         temp = ftrue(Qfull(:,k-1),time,time_end);
        F(sk,k) = ftrueSk(Q(:,k),time,time_end,sk');
        %         temp = ftrueSk(Q(:,k),time,time_end);
        %         F(sk,k) = temp(sk);
        F(skhat,k) = Uk(skhat,:)*pinv(Uk(sk,:))*F(sk,k);
    end
    
    time = time_end;
    
    %     sk = 1:m;
    %     Fk = Qfull(:, k-w+1:k);
    Fk = F(:, k-w+1:k);
    [Uk, Pk, ~] = adeim(Uk, Pk, sk, Fk(Pk,:), Fk(sk,:), r);
    %     [~,s,~] = svd(Fk,0);
    %     fprintf('d(Uk+1, Ubark+1) = %e\n', rho2/(min(diag(s))^2));
    
    
    if(k == wtotal)
        %% Local coherence
%         qnew = ftilde(Q(:,k),time,time_end+tstep,Uk,Pk); % Ck = (Pk'*Uk)\Pk'*Fk
%         F = ftrue(Uk*qnew,time,time_end+tstep);
        F = Qfull(:,k+1);
        R = -F + Uk*(Uk(Pk,:)\F(Pk));
        [~,idx] = sort(sum((R.^2),2),'descend');
        ressqr = R(idx).^2;
        figure(3);
        semilogy(ressqr, '-k');
        save('res10.mat','ressqr');
        save('restime10.mat','time');
    end
end
toc
%%
figure(4);
semilogy(1:wtotal-(winit),(errs'),'o-','Linewidth',1.5);
title("|q_{true}-UU^Tq_{true}|");
ylabel('error');