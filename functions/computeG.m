function [g] = computeG(S,psi)    
    S_H = ctranspose(S);
    SH_invPsi_S = S_H * inv(psi) * S;
    inv_SH_invPsi_S = inv(SH_invPsi_S);
    g = sqrt( diag(inv_SH_invPsi_S).*diag(SH_invPsi_S) );
end

