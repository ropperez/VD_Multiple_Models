function [z_series] = get_measurements_generation(X_truth,chol_R,H_gen,N_steps,N_iter)
% Measurements at k = 1
z_series = zeros(2,N_steps,N_iter);

for n1 = 1:N_iter
    for n2 = 1:N_steps
        % Measurement
        z_series(:,n2,n1) = H_gen*X_truth(:,n2)+chol_R*randn(2,1);
    end
end
end