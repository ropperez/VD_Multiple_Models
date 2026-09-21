function [logLr_k] = get_loglikelihood(z_k_inn,S_k)
m       = size(z_k_inn,1); % gaussian dimension
log_c   = (-m/2)*log(2*pi) + (-1/2)*log(det(S_k));
log_exp = -(1/2)*z_k_inn'/S_k*z_k_inn;
% log_exp = -(1/2)*z_k_inn'*pinv(S_k)*z_k_inn;
logLr_k = log_c + log_exp;
end