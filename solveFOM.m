function [Qh, Qv, time] = solveFOM(h, v, B, time, tstep, winit, limiter)
    [Np, K] = size(h);
    Qh = zeros(Np*K, winit);
    Qv = zeros(Np*K, winit);
    for i = 1:winit
        time_end = time + tstep;
        [h, v] = StateFixTS1D(h, v, time, time_end, B, limiter);
        time = time_end;
        Qh(:,i) = reshape(h, [Np*K, 1]);
        Qv(:,i) = reshape(v, [Np*K, 1]);
    end
end