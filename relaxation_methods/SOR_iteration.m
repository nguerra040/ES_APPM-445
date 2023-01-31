function [u] = SOR_iteration(u, N, eps, w)
% Author: Nicolas Guerra
% Date: January 29, 2023
%
% This function performs one iteration of Successive Over Relaxation for 
% the diff eq: -(epsilon*u_xx+u_yy)=0 with zero at the boundaries
% This method uses immediate replacement.
% INPUT: (matrix) u, (int) N, (double) eps, (double) w
% OUTPUT: (matrix) u
    for j = 2:(N+1)
        for k = 2:(N+1)
            Au = eps*u(j+1,k)+eps*u(j-1,k)-(2*eps+2)*u(j,k)+u(j,k+1)+u(j,k-1);
            u(j,k) = -w*Au/(-(2*eps+2)) + u(j,k);
        end
    end
end