function [E, iter] = kepler_solve(M, e, tol, E0)
% KEPLER_SOLVE Solves Kepler's equation M = E - e*sin(E) using Newton-Raphson.
%
% Inputs:
% M= Mean anomaly (rad)
% e=  Orbital eccentricity (-)
% tol= Convergence tolerance (rad)
% E0= Initial guess for eccentric anomaly (rad)
%
% Outputs:
% E = Converged eccentric anomaly (rad)
% iter= Number of iterations taken

max_iter = 50; 
E = E0;
iter = 0;

for k = 1:max_iter
    iter = k;
    f_val  = E - e * sin(E) - M;      
    fp_val = 1 - e * cos(E);          

    dE = -f_val / fp_val;             % Newton-Raphson update step
    E  = E + dE;

    % Step-size convergence criterion: |E_{k+1} - E_k| < tol
    if abs(dE) < tol
        return;
    end
end
end

%[appendix]{"version":"1.0"}
%---
