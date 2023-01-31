function [u_init] = starting_data_a(N)
% Author: Nicolas Guerra
% Date: January 29, 2023
%
% This function takes in a value N and returns a matrix of size (N+2)^2
% containing the initial data as specific in (a) of relaxation project
% INPUT: (int) N
% OUTPUT: (matrix of size (N+2)^2) u
    u_init = zeros(N+2,N+2);
    for j = 2:(N+1)
        for k = 2:(N+1)
            u_init(j,k) = (-1)^(j+k);
        end
    end
end