N = 16; % Number of interior points
eps = .1;
tol = 10^-7;

% (1) initial u 
u = zeros(N+2,N+2);
for j = 2:(N+1)
    for k = 2:(N+1)
        u(j,k) = (-1)^(j+k);
    end
end

% % Point Jacobi
% num_iter = 0;
% while max(max(abs(u))) > tol
%     u = point_jacobi_iteration(u, N, eps);
%     num_iter = num_iter + 1;
% end

% % Weighted Jacobi
% num_iter = 0;
% w = 2/3;
% while max(max(abs(u))) > tol
%     u = weighted_jacobi_iteration(u, N, eps, w);
%     num_iter = num_iter + 1;
% end

% % Gauss-Seidel
% num_iter = 0;
% while max(max(abs(u))) > tol
%     u = gauss_seidel_iteration(u, N, eps);
%     num_iter = num_iter + 1;
% end

% % Red-Black Gauss-Seidel
% num_iter = 0;
% while max(max(abs(u))) > tol
%     u = rb_gauss_seidel_iteration(u, N, eps);
%     num_iter = num_iter + 1;
% end

% % SOR
% num_iter = 0;
% w = 1.8;
% while max(max(abs(u))) > tol
%     u = SOR_iteration(u, N, eps, w);
%     num_iter = num_iter + 1;
% end

% % SSOR
% num_iter = 0; % remember, one iteration actually includes two sweeps
% w = 1;
% while max(max(abs(u))) > tol
%     u = SSOR_iteration(u, N, eps, w);
%     num_iter = num_iter + 1;
% end



% Kaczmarz
% Construct matrix for Kaczmarz for interior points
A = diag(-ones(N^2,1)*(2*eps+2))+...
    diag(ones(N^2-1,1)*eps,1)+...
    diag(ones(N^2-1,1)*eps,-1)+...
    diag(ones(N^2-N,1),N)+...
    diag(ones(N^2-N,1),-N);
for i = N:N:(N^2-N)
    A(i,i+1) = 0;
    A(i+1,i) = 0;
end
A = sparse(A);

% for i = 1:100000
%     u = kaczmarz_iteration(u, A);
% end
% disp(max(max(abs(u))))








