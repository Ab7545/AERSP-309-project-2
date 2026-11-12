
clear; clc; close all;

r1 = [3420.8; -9953.8; -657.2];   
v1 = [3.8415; 4.2736; -4.5964];  
mu = 3.986e5;                     

r1_mag = norm(r1);
v1_mag = norm(v1);


e_vec = (1/mu) * ((v1_mag^2 - mu/r1_mag)*r1 - dot(r1, v1)*v1); % Eccentricity vector
e_orbit = norm(e_vec);

fprintf('Project Orbit Eccentricity e = %.6f\n', e_orbit); %[output:4f4db3ba]

N_M      = 1000;
M_vals   = linspace(0, 2*pi*(1 - 1/N_M), N_M);
tol      = 1e-5;                               
max_iter = 50;                                 

e_cases     = [e_orbit, 0.99];
e_titles    = {sprintf('Project Orbit (e = %.4f)', e_orbit), ...
    'High-Eccentricity Orbit (e = 0.99)'};
guess_names = {'E_0 = M', 'E_0 = M + e sin(M)', 'E_0 = \pi'};

%%\ Newton-Raphson Across Both Eccentricities & All 3 Initial Guesses
figure('Name', 'Problem 1(d): Newton-Raphson Convergence', 'Position', [100 100 850 650]); %[output:0d46d2b8]

for idx_e = 1:2 %[output:group:5a23a7a5]
    ec = e_cases(idx_e);
    iters  = zeros(3, N_M);
    failed = false(3, N_M);

    for j = 1:N_M
        M = M_vals(j);
        E0_list = [M, M + ec*sin(M), pi];

        for g = 1:3
            [E_sol, it] = kepler_solve(M, ec, tol, E0_list(g));
            iters(g, j) = it;

            % Checking  if the  solver failed to meet tolerance within 50 iterations
            if (it >= max_iter && abs(E_sol - ec*sin(E_sol) - M) >= tol) || isnan(E_sol) || isinf(E_sol)
                failed(g, j) = true;
            end
        end
    end

    %%\ 4. Plot Iterations vs. M (Logarithmic Vertical Axis)
    subplot(2, 1, idx_e); %[output:0d46d2b8]
    semilogy(M_vals, iters(1, :), 'b-',  'LineWidth', 1.2); hold on;
    semilogy(M_vals, iters(2, :), 'r--', 'LineWidth', 1.2);
    semilogy(M_vals, iters(3, :), 'g-.', 'LineWidth', 1.2);
    grid on;
    xlim([0, 2*pi]);
    xticks([0, pi/2, pi, 3*pi/2, 2*pi]);
    xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    xlabel('Mean Anomaly M (rad)');
    ylabel('Number of Iterations (-)');
    title(['Newton-Raphson Iterations vs. Mean Anomaly — ' e_titles{idx_e}]);
    legend(guess_names, 'Location', 'north');

    %%\ 5. Print Summary Statistics Table for Report
    fprintf('\n====================================================================\n'); %[output:1d8738b8] %[output:0b6a3e8c]
    fprintf('Case: %s (tol = %.0e, max_iter = %d)\n', e_titles{idx_e}, tol, max_iter); %[output:1469d0ca] %[output:1eaeeceb]
    fprintf('====================================================================\n'); %[output:9da49482] %[output:742edce5]
    fprintf('%-22s | %-12s | %-12s | %-10s\n', 'Initial Guess', 'Mean Iters', 'Max Iters', 'Failures'); %[output:1bfd2047] %[output:58638fb7]
    fprintf('--------------------------------------------------------------------\n'); %[output:2b28a964] %[output:0c1204e2]
    for g = 1:3
        fprintf('%-22s | %12.2f | %12d | %10d\n', ... %[output:5217d5a2] %[output:7bda2eec]
            guess_names{g}, mean(iters(g, :)), max(iters(g, :)), sum(failed(g, :))); %[output:5217d5a2] %[output:7bda2eec]
    end
end 
