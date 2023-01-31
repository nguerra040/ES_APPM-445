function [u] = gauss_seidel_iteration(u, N, eps)
% Author: Nicolas Guerra
% Date: January 29, 2023
%
% This function performs one iteration of Gauss-Seidel for 
% the diff eq: -(epsilon*u_xx+u_yy)=0 with zero at the boundaries
% INPUT: (matrix) u, (int) N, (double) eps
% OUTPUT: (matrix) u
    for j = 2:(N+1)
        for k = 2:(N+1)
            % Doesn't include u you are currently updating in equation
            Au = eps*u(j+1,k)+eps*u(j-1,k)+u(j,k+1)+u(j,k-1);
            u(j,k) = (-Au)/(-(2*eps+2));
        end
    end
end