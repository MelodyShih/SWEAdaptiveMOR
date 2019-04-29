function [P] = qdeim(U)
    [~, ~, P] = qr(U', 'vector');
    P = P(1:size(U,2));
end

