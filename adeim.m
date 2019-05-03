function [Unew,Pnew] = adeim(U,P,S,Fp,Fs,r)
    C = U(P,:)\Fp; %% ftilde
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
    [Q,rtest] = qr(U); %Orthogonalization of U
    Unew = Q(:,1:size(U,2));
    Pnew = qdeim(Unew);
    % Test whether the new basis can represent the state 
    %Cnew = Unew(P,:)\Fp; %% ftilde
    %Rnew = Unew(S,:)*Cnew-Fs;
    %fprintf('original max error: %e, new basis: %e\n', ...
    %   max(max(R)), max(max(Rnew)));
end