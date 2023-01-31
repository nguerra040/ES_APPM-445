function [new_u] = point_jacobi_iteration(u, N, eps)
% Author: Nicolas Guerra
% Date: January 6, 2023
%
% This function performs one iteration of point Jacobi for 
% the diff eq: -(epsilon*u_xx+u_yy)=0 with zero at the boundaries
% INPUT: (matrix) u, (int) N, (double) eps
% OUTPUT: (matrix) u
    new_u = u;
    for j = 2:(N+1)
        for k = 2:(N+1)
            Au = eps*u(j+1,k)+eps*u(j-1,k)-(2*eps+2)*u(j,k)+u(j,k+1)+u(j,k-1);
            new_u(j,k) = (-(2*eps+2)*u(j,k) - Au)/(-(2*eps+2));
        end
    end
end