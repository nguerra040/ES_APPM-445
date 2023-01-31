function [optimal_w] = get_optimal_w(method, u_init, N, eps)

if method == "Weighted Jacobi"
    w = 0.5:.01:1;
    errors = zeros(length(w),1);
    for i = 1:length(w) % cycle through each w
        u = u_init;
        for j = 1:100 % relax for 100 iterations
            u = weighted_jacobi_iteration(u, N, eps, w(i));
        end
        errors(i) = max(max(abs(u)));
    end
    % Return optimal w
    [~,index] = min(errors);
    optimal_w = w(index);


elseif method == "SOR"
    w = 1.5:.01:1.99;
    errors = zeros(length(w),1);
    for i = 1:length(w) % cycle through each w
        u = u_init;
        for j = 1:100 % relax for 100 iterations
            u = SOR_iteration(u, N, eps, w(i));
        end
        errors(i) = max(max(abs(u)));
    end
    % Return optimal w
    [~,index] = min(errors);
    optimal_w = w(index);

elseif method == "SSOR"
    w = 1.5:.01:1.99;
    errors = zeros(length(w),1);
    for i = 1:length(w) % cycle through each w
        u = u_init;
        for j = 1:100 % relax for 100 iterations
            u = SSOR_iteration(u, N, eps, w(i));
        end
        errors(i) = max(max(abs(u)));
    end
    % Return optimal w
    [~,index] = min(errors);
    optimal_w = w(index);

end

end