function [U] = computeUnfolding(S,Gamma)
    % U = (S^H * gamma^-1 * S)^-1 * S^H *gamma^-1
    U = (ctranspose(S)/Gamma*S)\ctranspose(S)/Gamma;
end

