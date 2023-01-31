% Author: Nicolas Guerra
% Date: January 29, 2023
%
% This program iterates Weighted Jacobi, SOR, and SSOR 100 times to solve
% the PDE of interest using varying values of omega. After cycling through
% the various values of omega, a plot is graphed of the max residual as a
% function of omega. This shows the sensitivity of convergence as a
% function of omega.

N = 128; % Number of interior points
eps = 1;
tol = 10^-7;

% Weighted Jacobi
w_wj = 0.91:0.01:.95;
max_residual_wj = zeros(size(w_wj));
for i = 1:length(w_wj)
    % (a) initial u 
    u = starting_data_a(N);

    % Weighted Jacobi 100 iterations
    for j = 1:100
        u = weighted_jacobi_iteration(u, N, eps, w_wj(i));
    end
    % Save number of iterations
    max_residual_wj(i) = max(max(abs(u)));
end
figure(1)
plot(w_wj,max_residual_wj)
xlabel('\omega')
ylabel('max|u_{ij}|')

% SOR
w_SOR = 1.71:0.01:1.99;
max_residual_SOR = zeros(size(w_SOR));
for i = 1:length(w_SOR)
    u = starting_data_c(N);
        % SOR 100 iterations
    for j = 1:100
        u = SOR_iteration(u, N, eps, w_SOR(i));
    end
    % Save number of iterations
    max_residual_SOR(i) = max(max(abs(u)));
end
figure(2)
plot(w_SOR,max_residual_SOR)
xlabel('\omega')
ylabel('max|u_{ij}|')

% SSOR
w_SSOR = 1.71:0.01:1.99;
max_residual_SSOR = zeros(size(w_SSOR));
for i = 1:length(w_SSOR)
    u = starting_data_c(N);
        % SOR 100 iterations
    for j = 1:100
        u = SSOR_iteration(u, N, eps, w_SSOR(i));
    end
    % Save number of iterations
    max_residual_SSOR(i) = max(max(abs(u)));
end
figure(3)
plot(w_SSOR,max_residual_SSOR)
xlabel('\omega')
ylabel('max|u_{ij}|')