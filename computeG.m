function [g] = computeG(S,Gamma)    
    g = diag(sqrt((ctranspose(S)/Gamma*S)\(ctranspose(S)/Gamma*S)));
end

