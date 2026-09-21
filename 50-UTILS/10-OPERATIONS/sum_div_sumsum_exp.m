function [y] = sum_div_sumsum_exp(z_k_inn_ji,S_k_ji,w_ji,N_r)
% Pre-allocation memory
log_Lr_ji = zeros(size(z_k_inn_ji));

% Get likelihood
for i = 1:N_r
    for j = 1:N_r
        log_Lr_ji(j,i) = get_loglikelihood(z_k_inn_ji{j,i},S_k_ji{j,i});
    end
end

% Weighted log-Likelihood
log_w_Lr_ji = log(w_ji) + log_Lr_ji;

% Pre-allocation 
L = zeros(1,N_r);

% Output in natural units
y = zeros(1,N_r);
for i = 1:N_r
    L(i) = max(log_w_Lr_ji(:,i));
    y(i) = y(i) + sum(exp(log_w_Lr_ji(:,i)-L(i)))/sum(sum(exp(log_w_Lr_ji-L(i)),1));
end

end