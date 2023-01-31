% Author: Nicolas Guerra
% Date: January 29, 2023
%
% This program goes through every combination of N and epsilon for each
% starting data and solves the PDE of interest using six different methods. 
% It then stores the number of iterations it takes to reach the
% desired tolerance

N = [16, 32, 64, 128]; % Number of interior points
eps = [1, .1, .01, .001];
tol = 10^-7;

%% Starting Data (a)
sz = [112 5];
var_types = ["string","double","double","double","double"];
var_names = ["Method","N","Epsilon","Omega","Number of Iterations"];
output_a = table('Size',sz,'VariableTypes',var_types,'VariableNames',var_names);
fill_row = 1;
for i = 1:length(N)
    % (a) initial u 
    u_init = starting_data_a(N(i));

    for j = 1:length(eps)
        % Point Jacobi Loop
        [~, num_iter] = point_jacobi_loop(u_init, tol, N(i), eps(j));
        % Save number of iterations
        output_a(fill_row,:)={'Point Jacobi',N(i),eps(j),NaN,num_iter};
        fill_row = fill_row + 1;


        % Weighted Jacobi Loop
        w = get_optimal_w('Weighted Jacobi', u_init, N(i), eps(j));
        [~, num_iter] = weighted_jacobi_loop(u_init, tol, N(i), eps(j), w);
        % Save number of iterations
        output_a(fill_row,:)={'Weighted Jacobi',N(i),eps(j),w,num_iter};
        fill_row = fill_row + 1;

        % Gauss-Seidel Loop
        [~, num_iter] = gauss_seidel_loop(u_init, tol, N(i), eps(j));
        % Save number of iterations
        output_a(fill_row,:)={'Gauss-Seidel',N(i),eps(j),NaN,num_iter};
        fill_row = fill_row + 1;

        % Red-Black Gauss-Seidel Loop
        [~, num_iter] = rb_gauss_seidel_loop(u_init, tol, N(i), eps(j));
        % Save number of iterations
        output_a(fill_row,:)={'Red-Black Gauss-Seidel',N(i),eps(j),NaN,num_iter};
        fill_row = fill_row + 1;

        % SOR Loop
        w = get_optimal_w('SOR', u_init, N(i), eps(j));
        [~, num_iter] = SOR_loop(u_init, tol, N(i), eps(j), w);
        % Save number of iterations
        output_a(fill_row,:)={'SOR',N(i),eps(j),w,num_iter};
        fill_row = fill_row + 1;

        % SSOR Loop
        w = get_optimal_w('SSOR', u_init, N(i), eps(j));
        [~, num_iter] = SSOR_loop(u_init, tol, N(i), eps(j), w);
        % Save number of iterations
        output_a(fill_row,:)={'SSOR',N(i),eps(j),w,num_iter};
        fill_row = fill_row + 1;

        % Kaczmarz Loop
        [~, num_iter] = kaczmarz_loop(u_init, tol, N(i), eps(j));
        % Save number of iterations
        output_a(fill_row,:)={'Kaczmarz',N(i),eps(j),NaN,num_iter};
        fill_row = fill_row + 1;
    end
end

%% Starting Data (b)
output_b = table('Size',sz,'VariableTypes',var_types,'VariableNames',var_names);
fill_row = 1;
for i = 1:length(N)
    % (b) initial u 
    u_init = starting_data_b(N(i));

    for j = 1:length(eps)
        % Point Jacobi Loop
        [~, num_iter] = point_jacobi_loop(u_init, tol, N(i), eps(j));
        % Save number of iterations
        output_b(fill_row,:)={'Point Jacobi',N(i),eps(j),NaN,num_iter};
        fill_row = fill_row + 1;


        % Weighted Jacobi Loop
        w = get_optimal_w('Weighted Jacobi', u_init, N(i), eps(j));
        [~, num_iter] = weighted_jacobi_loop(u_init, tol, N(i), eps(j), w);
        % Save number of iterations
        output_b(fill_row,:)={'Weighted Jacobi',N(i),eps(j),w,num_iter};
        fill_row = fill_row + 1;

        % Gauss-Seidel Loop
        [~, num_iter] = gauss_seidel_loop(u_init, tol, N(i), eps(j));
        % Save number of iterations
        output_b(fill_row,:)={'Gauss-Seidel',N(i),eps(j),NaN,num_iter};
        fill_row = fill_row + 1;

        % Red-Black Gauss-Seidel Loop
        [~, num_iter] = rb_gauss_seidel_loop(u_init, tol, N(i), eps(j));
        % Save number of iterations
        output_b(fill_row,:)={'Red-Black Gauss-Seidel',N(i),eps(j),NaN,num_iter};
        fill_row = fill_row + 1;

        % SOR Loop
        w = get_optimal_w('SOR', u_init, N(i), eps(j));
        [~, num_iter] = SOR_loop(u_init, tol, N(i), eps(j), w);
        % Save number of iterations
        output_b(fill_row,:)={'SOR',N(i),eps(j),w,num_iter};
        fill_row = fill_row + 1;

        % SSOR Loop
        w = get_optimal_w('SSOR', u_init, N(i), eps(j));
        [~, num_iter] = SSOR_loop(u_init, tol, N(i), eps(j), w);
        % Save number of iterations
        output_b(fill_row,:)={'SSOR',N(i),eps(j),w,num_iter};
        fill_row = fill_row + 1;

        % Kaczmarz Loop
        [~, num_iter] = kaczmarz_loop(u_init, tol, N(i), eps(j));
        % Save number of iterations
        output_b(fill_row,:)={'Kaczmarz',N(i),eps(j),NaN,num_iter};
        fill_row = fill_row + 1;
    end
end
 
%% Starting Data (c)
output_c = table('Size',sz,'VariableTypes',var_types,'VariableNames',var_names);
fill_row = 1;
for i = 1:length(N)
    % (c) initial u 
    u_init = starting_data_c(N(i));

    for j = 1:length(eps)
        % Point Jacobi Loop
        [~, num_iter] = point_jacobi_loop(u_init, tol, N(i), eps(j));
        % Save number of iterations
        output_c(fill_row,:)={'Point Jacobi',N(i),eps(j),NaN,num_iter};
        fill_row = fill_row + 1;


        % Weighted Jacobi Loop
        w = get_optimal_w('Weighted Jacobi', u_init, N(i), eps(j));
        [~, num_iter] = weighted_jacobi_loop(u_init, tol, N(i), eps(j), w);
        % Save number of iterations
        output_c(fill_row,:)={'Weighted Jacobi',N(i),eps(j),w,num_iter};
        fill_row = fill_row + 1;

        % Gauss-Seidel Loop
        [~, num_iter] = gauss_seidel_loop(u_init, tol, N(i), eps(j));
        % Save number of iterations
        output_c(fill_row,:)={'Gauss-Seidel',N(i),eps(j),NaN,num_iter};
        fill_row = fill_row + 1;

        % Red-Black Gauss-Seidel Loop
        [~, num_iter] = rb_gauss_seidel_loop(u_init, tol, N(i), eps(j));
        % Save number of iterations
        output_c(fill_row,:)={'Red-Black Gauss-Seidel',N(i),eps(j),NaN,num_iter};
        fill_row = fill_row + 1;

        % SOR Loop
        w = get_optimal_w('SOR', u_init, N(i), eps(j));
        [~, num_iter] = SOR_loop(u_init, tol, N(i), eps(j), w);
        % Save number of iterations
        output_c(fill_row,:)={'SOR',N(i),eps(j),w,num_iter};
        fill_row = fill_row + 1;

        % SSOR Loop
        w = get_optimal_w('SSOR', u_init, N(i), eps(j));
        [~, num_iter] = SSOR_loop(u_init, tol, N(i), eps(j), w);
        % Save number of iterations
        output_c(fill_row,:)={'SSOR',N(i),eps(j),w,num_iter};
        fill_row = fill_row + 1;

        % Kaczmarz Loop
        [~, num_iter] = kaczmarz_loop(u_init, tol, N(i), eps(j));
        % Save number of iterations
        output_c(fill_row,:)={'Kaczmarz',N(i),eps(j),NaN,num_iter};
        fill_row = fill_row + 1;
    end
end