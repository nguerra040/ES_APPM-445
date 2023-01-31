function [u] = rb_gauss_seidel_iteration(u, N, eps)
% Author: Nicolas Guerra
% Date: January 29, 2023
%
% This function performs one iteration of Red-Black Gauss-Seidel for 
% the diff eq: -(epsilon*u_xx+u_yy)=0 with zero at the boundaries
% INPUT: (matrix) u, (int) N, (double) eps
% OUTPUT: (matrix) u
for j = 2:(N+1)
    for k = 2:2:(N+1)
        zigzag_k = k+mod(j,2);
        % Black
        % Don't include u you are currently updating in equation
        Au = eps*u(j+1,zigzag_k)+eps*u(j-1,zigzag_k)+u(j,zigzag_k+1)+u(j,zigzag_k-1);
        u(j,zigzag_k) = (-Au)/(-(2*eps+2));
    end
end
for j = 2:(N+1)
    for k = 2:2:(N+1)
        zigzag_k = k+mod(j+1,2);
        % Red
        % Don't include u you are currently updating in equation
        Au = eps*u(j+1,zigzag_k)+eps*u(j-1,zigzag_k)+u(j,zigzag_k+1)+u(j,zigzag_k-1);
        u(j,zigzag_k) = (-Au)/(-(2*eps+2));
    end
end