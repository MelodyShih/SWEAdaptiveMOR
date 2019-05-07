function [Krd,crd,crdm] = ReducePoint(Sk)
% Build the map from the reduced basis locations to the full space
% Input: Sk -- a 2Np*K row vector gives the location of the points
% Output: Krd -- the number of columns we use 
%         crd -- the column number we use for vector form
%         crdm -- the column number we use for matrix form

Globals1D;

totalpts = Np*K; % total number of points
crd = ceil((Sk.*(Sk<=totalpts)+(Sk-totalpts).*(Sk>totalpts))/Np);
crd = unique([crd+1,crd,crd-1]);
crd(crd<=1)=[];
crd(crd>=K)=[];
crdm=[1,crd,K]; % columns for matrix form
Krd = length(crdm);
crd = zeros(1,2*Krd); % vector form
crd(1:2:end)=crdm*2-1;
crd(2:2:end)=crdm*2;
end

