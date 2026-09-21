function [y] = div_sum_exp(z_k_inn_i,S_k_i,w_i,N_r)
% Memory preallocation
log_Lr_k = zeros(N_r,1); 

% Get likelihood
for i = 1:N_r
    log_Lr_k(i) = get_loglikelihood(z_k_inn_i(:,i),S_k_i(:,:,i));
end

% Weigthed log-Likelihood
log_w_Lr_i  = log(w_i) + log_Lr_k;

% Max weigthed likelihood (constant)
L           = max(log_w_Lr_i);

% Output in natural units
y           = exp(log_w_Lr_i-L)/sum(exp(log_w_Lr_i-L));
end
