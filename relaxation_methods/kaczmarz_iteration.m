function [u_new] = kaczmarz_iteration(u,A)
% Author: Nicolas Guerra
% Date: January 29, 2023
%
% This function performs one iteration of Kaczmarz for 
% the diff eq: -(epsilon*u_xx+u_yy)=0 with zero at the boundaries
% INPUT: (matrix) u, (matrix) A
% OUTPUT: (matrix) u_new
    A = sparse(A);
    u=u(2:end-1,2:end-1);
    N = size(u,1); % number of interior points
    x = reshape(u, N^2 ,1);
    for i = 1:N^2
        r_i = -A(i,:)*x; % residual of ith row
        delta_i = r_i/norm(A(i,:))^2;
        x = x + delta_i*A(i,:)';
    end
    reshaped_x = reshape(x, N, N);
    u_new = zeros(size(reshaped_x)+2);
    u_new(2:end-1,2:end-1)=reshaped_x;
end