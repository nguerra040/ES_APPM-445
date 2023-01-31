function [u, num_iter] = point_jacobi_loop(u_init, tol, N, eps)
% Author: Nicolas Guerra
% Date: January 29, 2023
%
% This function performs point Jacobi iterations until the desired
% tolerance is reached.
% INPUT: (matrix) u_init, (double) tol, (int) N, (double) eps
% OUTPUT: (matrix) u, (int) num_iter
    num_iter = 0;
    u = u_init;
    while max(max(abs(u))) > tol
        u = point_jacobi_iteration(u, N, eps);
        num_iter = num_iter + 1;
    end
end