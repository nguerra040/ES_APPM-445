% Author: Nicolas Guerra
% Date: January 29, 2023
%
% This program runs 10 iterations of the Kaczmarz method cyclying through
% each value of N and epsilon. The max residual after 10 iterations is then
% stored.

N = [16,32,64,128];
eps = [1,.1,.01,.001];

% Data (a)
max_residuals_a = [];
for i = 1:length(N)
    for j = 1:length(eps)
        u = starting_data_a(N(i));
        % Construct matrix for Kaczmarz for interior points
        A = diag(-ones(N(i)^2,1)*(2*eps(j)+2))+...
            diag(ones(N(i)^2-1,1)*eps(j),1)+...
            diag(ones(N(i)^2-1,1)*eps(j),-1)+...
            diag(ones(N(i)^2-N(i),1),N(i))+...
            diag(ones(N(i)^2-N(i),1),-N(i));
        for ii = N(i):N(i):(N(i)^2-N(i))
            A(ii,ii+1) = 0;
            A(ii+1,ii) = 0;
        end
        % Since A is mainly zeros, store as sparse matrix
        A = sparse(A);
        % Kaczmarz method
        for iter = 1:10
            u = kaczmarz_iteration(u, A);
        end
        % Store max residual
        max_residuals_a = [max_residuals_a; max(max(abs(u)))];
    end
end

% Data (b)
max_residuals_b = [];
for i = 1:length(N)
    for j = 1:length(eps)
        u = starting_data_b(N(i));
        % Construct matrix for Kaczmarz for interior points
        A = diag(-ones(N(i)^2,1)*(2*eps(j)+2))+...
            diag(ones(N(i)^2-1,1)*eps(j),1)+...
            diag(ones(N(i)^2-1,1)*eps(j),-1)+...
            diag(ones(N(i)^2-N(i),1),N(i))+...
            diag(ones(N(i)^2-N(i),1),-N(i));
        for ii = N(i):N(i):(N(i)^2-N(i))
            A(ii,ii+1) = 0;
            A(ii+1,ii) = 0;
        end
        % Since A is mainly zeros, store as sparse matrix
        A = sparse(A);
        % Kaczmarz method
        for iter = 1:10
            u = kaczmarz_iteration(u, A);
        end
        % Store max residual
        max_residuals_b = [max_residuals_b; max(max(abs(u)))];
    end
end

% Data (c)
max_residuals_c = [];
for i = 1:length(N)
    for j = 1:length(eps)
        u = starting_data_c(N(i));
        % Construct matrix for Kaczmarz for interior points
        A = diag(-ones(N(i)^2,1)*(2*eps(j)+2))+...
            diag(ones(N(i)^2-1,1)*eps(j),1)+...
            diag(ones(N(i)^2-1,1)*eps(j),-1)+...
            diag(ones(N(i)^2-N(i),1),N(i))+...
            diag(ones(N(i)^2-N(i),1),-N(i));
        for ii = N(i):N(i):(N(i)^2-N(i))
            A(ii,ii+1) = 0;
            A(ii+1,ii) = 0;
        end
        % Since A is mainly zeros, store as sparse matrix
        A = sparse(A);
        % Kaczmarz method
        for iter = 1:10
            u = kaczmarz_iteration(u, A);
        end
        % Store max residual
        max_residuals_c = [max_residuals_c; max(max(abs(u)))];
    end
end