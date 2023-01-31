% Author: Nicolas Guerra
% Date: January 29, 2023
%
% This program cycles through all the values of N and epsilon for a given
% the PDE of interest with the given initial data which can be set in line
% 15. The PDE is solved with SSOR with omega=1 to get symmetrized
% Gauss-Seidel.

% Values of N and epsilon to cycle through
N = [16,32,64,128];
eps = [1,.1,.01,.001];

num_iters = [];
for i = 1:length(N)
    for j = 1:length(eps)
        % Set initial data
        u = starting_data_c(N(i));
        % remember, one iteration actually includes two sweeps
        num_iter = 0; 
        % Set omega = 1 to get symmetrized Gauss-Seidel
        w = 1;
        % Solve the PDE using SSOR
        while max(max(abs(u))) > tol
            u = SSOR_iteration(u, N(i), eps(j), w);
            num_iter = num_iter + 1;
        end
        % Store the number of iterations it took
        num_iters = [num_iters; num_iter];
    end
end
