function [u, num_iter] = weighted_jacobi_loop(u_init, tol, N, eps, w)
% Author: Nicolas Guerra
% Date: January 29, 2023
%
% This function performs weighted Jacobi iterations until the desired
% tolerance is reached.
% INPUT: (matrix) u_init, (double) tol, (int) N, (double) eps, (double) w
% OUTPUT: (matrix) u, (int) num_iter
    num_iter = 0;
    u = u_init;
    while max(max(abs(u))) > tol
        u = weighted_jacobi_iteration(u, N, eps, w);
        num_iter = num_iter + 1;
    end
end