function [dim_c] = get_dim_c(r_k)

enum_dyn_model = enum_dynamic_model();

if r_k == enum_dyn_model.ca 
    dim_c = 3;
elseif r_k == enum_dyn_model.cv || r_k == enum_dyn_model.ct_l || r_k == enum_dyn_model.ct_r
    dim_c = 2;
end

end