function [dim_r] = get_dim_r(r_k)
% Enum and constants
% enum_dyn_model  = get_dynamic_model();

% switch r_k
%     case enum_dyn_model.cv
%         dim_r = 4;
%     case enum_dyn_model.ca
%         dim_r = 6;
%     case {enum_dyn_model.ct_r, enum_dyn_model.ct_l}
%         dim_r = 5;
% end

% Dimension of the dynamic models cv; ct_r; ct_l; ca;
dim_r = [4 5 5 4]; 

dim_r = dim_r(r_k);

end