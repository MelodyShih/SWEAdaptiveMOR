clear;
close all;

% Test the gradient
Globals1D;

% Final time
FinalTime = 0;

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


% initial condition
h = ones(size(x));

% setup bathymetry (mu)
a0 = (a+b)/2-10;
b0 = (a+b)/2+10;
mu = 0.15;
p = b0 - a0;
B = -h+mu*(1 + cos(2*pi/p*(x - (a0+b0)/2))).*(x>a0 & x<b0);

v = zeros(Np,K);

% fix time step
CFL=0.1; g=9.8;
mindx = min(abs(x(2,:)-x(1,:)));
tstep = CFL*min(min(mindx./(abs(v./h)+sqrt(g*h))));

% % solve the forward equation
% [lh, lv,lu, ltime, Nstep] = State1D(h, v, FinalTime, B);

% fix time step
% tstep = 1; % time step
time =0;

figure(1);
plot(x(1,:),h(1,:)+B(1,:),'b','LineWidth',2);
hold on;
plot(x(1,:),B(1,:),'k','LineWidth',2);
title(['t=',num2str(time)],'Interpreter','Latex','FontSize',15);
xlabel('$x$', 'Interpreter','Latex','FontSize',15);
legend({'water height','bathymetry'}, 'Interpreter', 'latex','FontSize',15);
ylim([-1, 0.5]);
hold off;

figure(2);
plot(x(1,:),v(1,:),'b','LineWidth',2);
xlabel('$x$', 'Interpreter','Latex','FontSize',15);
title(['t=',num2str(time)],'Interpreter','Latex','FontSize',15);
legend({'momentum'}, 'Interpreter', 'latex','FontSize',15);
ylim([-0.5,0.5]);

while time<FinalTime
    if FinalTime-time<tstep
        time_end = FinalTime;
    else
        time_end = time+tstep;
    end
    [h, v] = StateFixTS1D(h, v, time, time_end, B,limiter);
    time = time_end;
    
    figure(1);
    plot(x(1,:),h(1,:)+B(1,:),'b','LineWidth',2);
    hold on;
    plot(x(1,:),B(1,:),'k','LineWidth',2);
    title(['t=',num2str(time)],'Interpreter','Latex','FontSize',15);
    xlabel('$x$', 'Interpreter','Latex','FontSize',15);
    legend({'water height','bathymetry'}, 'Interpreter', 'latex','FontSize',15);
    ylim([-1, 0.5]);
    hold off;
    
    figure(2);
    plot(x(1,:),v(1,:),'b','LineWidth',2);
    xlabel('$x$', 'Interpreter','Latex','FontSize',15);
    title(['t=',num2str(time)],'Interpreter','Latex','FontSize',15);
    legend({'momentum'}, 'Interpreter', 'latex','FontSize',15);
    ylim([-0.5,0.5]);
    
end


hinit = h;
vinit = v;
save('data/h0.mat','hinit');
save('data/v0.mat','vinit');
save('data/time0.mat','time');
