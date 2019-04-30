function [Unew,Pnew] = adeim(U,P,S,Fp,Fs,r)
    C = U(P,:)\Fp;
    R = U(S,:)*C-Fs;
    [~,Sv,Sr] = svd(R,0);
    Sv = diag(Sv);
    Ctpinv = pinv(C');
    r = min([r,length(Sv)]);
    for i = 1:r
        alpha = -R*Sr(:,i);
        beta = Ctpinv*Sr(:,i);
        U(S,:) = U(S,:)+alpha*beta';
    end
    [Q,R] = qr(U); %Orthogonalization of U
    Unew = Q; 
    Pnew = qdeim(Unew); 
end

