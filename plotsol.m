function plotsol(x, q, qtrue, time)
    global Np K B
    h = reshape(q(1:Np*K), [Np, K]);
    v = reshape(q(Np*K+1:end), [Np, K]); % evaluate at all points first
    
    htrue = reshape(qtrue(1:Np*K), [Np, K]);
    vtrue = reshape(qtrue(Np*K+1:end), [Np, K]); % evaluate at all points first
    
    figure(1);
    plot(x(1,:),h(1,:)+B(1,:),'b','LineWidth',2);
    hold on;
    plot(x(1,:),htrue(1,:)+B(1,:)-1,'r','LineWidth',2);
    title('$h(x,t)+B(x)$', 'Interpreter', 'latex', 'Fontsize',20);
    xlabel('$x$', 'Interpreter', 'latex', 'Fontsize',20)
    legend({'approximation', 'true'}, 'Interpreter', 'latex', 'Fontsize',20);
    hold off;
    
    
    figure(2);
    plot(x(1,:),v(1,:),'b','LineWidth',2);
    hold on;
    plot(x(1,:),vtrue(1,:)-1,'r','LineWidth',2);
    title('$v(x,t)$', 'Interpreter', 'latex', 'Fontsize',20);
    legend({'approximation', 'true'}, 'Interpreter', 'latex', 'Fontsize',20);
    hold off;
end

