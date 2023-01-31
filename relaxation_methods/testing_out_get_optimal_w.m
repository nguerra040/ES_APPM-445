N = 12;
eps = 1;
u_init = starting_data_a(N);
method='Weighted Jacobi';
w = get_optimal_w(method, u_init, N, eps)