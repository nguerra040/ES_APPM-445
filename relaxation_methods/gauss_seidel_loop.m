function [u, num_iter] = gauss_seidel_loop(u_init, tol, N, eps)
% Author: Nicolas Guerra
% Date: January 29, 2023
%
% This function performs Gauss-Seidel iterations until the desired
% tolerance is reached.
% INPUT: (matrix) u_init, (double) tol, (int) N, (double) eps
% OUTPUT: (matrix) u, (int) num_iter
    num_iter = 0;
    u = u_init;
    while max(max(abs(u))) > tol
        u = gauss_seidel_iteration(u, N, eps);
        num_iter = num_iter + 1;
    end
end