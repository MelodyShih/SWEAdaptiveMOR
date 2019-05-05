function [Unew,Pnew,rhosquare] = adeim(U,P,S,Fp,Fs,r)
    C = U(P,:)\Fp; %% ftilde
    R = U(S,:)*C-Fs; 
    [~,Sv,Sr] = svd(R,0);
    Sv = diag(Sv);
    plot(Sv, '-o');
    rhosquare = sum(Sv(2:end).^2);
    fprintf('reduction factor rhok = %e\n', rhosquare);
    waitforbuttonpress;
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
    Cnew = Unew(P,:)\Fp; 
    Rnew = Unew(S,:)*Cnew-Fs; % should here be Cnew or C? 
                              % (compare to lemma 1, this should be C?)
    fprintf('original max error: %e, new basis: %e\n', ...
      norm(R, 'fro'), norm(Rnew, 'fro'));
end