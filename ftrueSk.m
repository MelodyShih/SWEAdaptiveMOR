function [fQ] = ftrueSk(Q,time,time_end,Sk)
    global limiter Np K B
    
    h = reshape(Q(1:Np*K), [Np, K]);
    v = reshape(Q(Np*K+1:end), [Np, K]); % evaluate at all points first
    
%     [h, v] = StateFixTS1D(h, v, time, time_end, B, limiter);
%     fQ = [reshape(h, [Np*K, 1]);reshape(v,[Np*K, 1])]; 
    
    [Krd,crd,crdm] = ReducePoint(Sk);
    [fQ] = StateFixTSRD1D(h,v,time,time_end,B,limiter,Krd,crd,crdm,Sk);
end
