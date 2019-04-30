function [Qh, Qv, time] = solveFOM(hinit, vinit, time, tstep, winit)
    global limiter Np K B

    Qh = zeros(Np*K, winit);
    Qv = zeros(Np*K, winit);
    for i = 1:winit
        time_end = time + tstep;
        [h, v] = StateFixTS1D(hinit, vinit, time, time_end, B, limiter);
        time = time_end;
        Qh(:,i) = reshape(h, [Np*K, 1]);
        Qv(:,i) = reshape(v, [Np*K, 1]);
    end
end