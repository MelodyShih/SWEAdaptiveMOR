function plotsol(x, q, qtrue, time)
    global Np K B
    h = reshape(q(1:Np*K), [Np, K]);
    v = reshape(q(Np*K+1:end), [Np, K]); % evaluate at all points first
    
    htrue = reshape(qtrue(1:Np*K), [Np, K]);
    vtrue = reshape(qtrue(Np*K+1:end), [Np, K]); % evaluate at all points first
    
    figure(1);
    plot(x,h+B,'b','LineWidth',2);
    hold on;
    plot(x,htrue+B-1,'r','LineWidth',2);
    hold off;
    title(['t=',num2str(time)]);
    
    figure(2);
    plot(x,v,'b','LineWidth',2);
    hold on;
    plot(x,vtrue-1,'r','LineWidth',2);
    hold off;
    title(['t=',num2str(time)]);
end

