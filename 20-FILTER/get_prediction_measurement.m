function [z_k_hat,S_k] = get_prediction_measurement(x_k_hat,P_k_hat,H_k,R_k)
% Predicted measurement
z_k_hat = H_k*x_k_hat;
S_k     = H_k*P_k_hat*H_k'+R_k;
end