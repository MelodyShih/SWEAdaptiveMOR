%% Setup solver
% Test the gradient
Globals1D;

% Final time
FinalTime = 50;

% given parameter 
a=0; b=500; % space ind_k

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
mu_totry = [0.15];
CFL=0.2; g=9.8;
mindx = min(abs(x(2,:)-x(1,:)));
tstep = CFL*min(min(mindx./(abs(v./h)+sqrt(g*h))));
w = ceil(FinalTime/tstep); % store all the solutions
Qh = zeros(Np*K, w*length(mu_totry));
Qv = zeros(Np*K, w*length(mu_totry));


for i=1:length(mu_totry)
    fprintf('mu = %f starting...\n', mu_totry(i))
    % initial condition
    h = ones(size(x)); % size(x) = (#order of poly+1)x(#discretization pts)
    
    % setup bathymetry (mu)
    a0 = 220;
    b0 = 280;
    mu = mu_totry(i);
    p = b0 - a0;
    B = -h+mu*(1 + cos(2*pi/p*(x - (a0+b0)/2))).*(x>a0 & x<b0);
    v = zeros(Np,K);

    limiter = 1; % if use limiter

    %% Plot initial condition
    % fix time step 
    time =0;
    figure(1);
    plot(x,h+B,'b',x,B,'k','LineWidth',2);
    title(['t=',num2str(time)]);


    %% Start simulation
    total_step =1;
    while time<FinalTime
        if FinalTime-time<tstep
            time_end = FinalTime;
        else
            time_end = time+tstep;
        end
        [h, v] = StateFixTS1D(h, v, time, time_end, B,limiter);
        time = time_end;
        Qh(:,(i-1)*w+total_step) = reshape(h, [Np*K, 1]);
        Qv(:,(i-i)*w+total_step) = reshape(v, [Np*K, 1]);
        total_step =total_step +1;
%         figure(1);
%         plot(x,h+B,'b',x,B,'k','LineWidth',2);
%         title(['t=',num2str(time)]);
    end
end

%% local low-rankness
close all;
local_step = ceil(total_step/20);
local_start=1;
for i=1:length(mu_totry)
    s = svd(Qh(:, (i-1)*w+1:i*w) - mean(Qh(:, (i-1)*w+1:i*w)));
%     s_local = svd(Qh(:, (i-1)*w+1:(i-1)*w+local_step+1) - mean(Qh(:, (i-1)*w+1:(i-1)*w+local_step+1)));
    s_local = svd(Qh(:, (i-1)*w+local_start:(i-1)*w+local_step+local_start) - mean(Qh(:, (i-1)*w+local_start:(i-1)*w+local_step+local_start)));
    figure;
    semilogy(s/s(1), '-o')
    hold on;
    semilogy(s_local/s_local(1), '-*')
    title(['local low-rank structure, $\mu$=', num2str(mu_totry(i)), ', FinalTime=',num2str(FinalTime),', space dimension=',num2str(Np*K)], 'Interpreter', 'Latex')
    ylabel('normalized singular value', 'Interpreter', 'Latex')
    legend(['global, steps =',num2str(total_step)], ['local, steps =',num2str(local_step)], 'Interpreter', 'Latex')
    grid on;
end
