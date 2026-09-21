function [x_k, P_k] = get_update_stage(z_k_inn,S_k,x_k_hat,P_k_hat,H_k,R_k)
% Kalman Gain
K_k = P_k_hat*H_k'/S_k;

% State Update and associated covariance
x_k = x_k_hat+K_k*z_k_inn;
% P_k = P_k_hat-K_k*H_k*P_k_hat;   
P_k = (eye(size(P_k_hat))-K_k*H_k)*P_k_hat*(eye(size(P_k_hat))-K_k*H_k)'+K_k*R_k*K_k'; % Joseph Form
end