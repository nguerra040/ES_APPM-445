function [u, num_iter] = kaczmarz_loop(u_init, tol, N, eps)

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
    num_iter = 0; 
    u = u_init;
    while max(max(abs(u))) > tol
        u = kaczmarz_iteration(u, A);
        num_iter = num_iter + 1;
    end
end