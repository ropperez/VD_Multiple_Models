function [x_k_i,P_k_i] = get_init_state_cov(x_k,P_k,r_k,ENUM)
% Enum
enum_dyn_model = ENUM.dynamic_model;

if r_k == enum_dyn_model.cv
    x_k_i = x_k(1:4);
    P_k_i = P_k(1:4,1:4);
elseif r_k == enum_dyn_model.ct_l || r_k == enum_dyn_model.ct_r
    x_k_i = x_k;
    P_k_i = P_k;
end

end