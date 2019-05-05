function [Qh, Qv, time] = solveFOM(hinit, vinit, time, tstep, winit)
    global limiter Np K B

    Qh = zeros(Np*K, winit);
    Qv = zeros(Np*K, winit);
    
    h = hinit;
    v = vinit;
    for i = 1:winit
        time_end = time + tstep;
        [h, v] = StateFixTS1D(h, v, time, time_end, B, limiter);
        time = time_end;
        Qh(:,i) = reshape(h, [Np*K, 1]);
        Qv(:,i) = reshape(v, [Np*K, 1]);
%         if(i == 500)
%             hinit = h;
%             vinit = v;
%             save('h500.mat','hinit');
%             save('v500.mat','vinit');
%             save('time.mat','time');
%         end
    end
    
end