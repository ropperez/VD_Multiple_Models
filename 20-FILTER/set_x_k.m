function [x_k_dim] = set_x_k(x_k,r_k,dim)
enum_dyn_model = enum_dynamic_model();

x_k_dim = zeros(1,dim);

if dim == 7 || dim == 6
    if r_k == enum_dyn_model.ca
        x_k_dim(1:6) = x_k;
    elseif r_k == enum_dyn_model.cv || r_k == enum_dyn_model.ct_l || r_k == enum_dyn_model.ct_r
        x_k_dim(1:2) = x_k(1:2);
        x_k_dim(4:5) = x_k(3:4);
    end
elseif dim == 5
    if r_k == enum_dyn_model.cv
        x_k_dim(1:4) = x_k(1:4);
    elseif r_k == enum_dyn_model.ct_l || r_k == enum_dyn_model.ct_r
        x_k_dim = x_k;
    end    
end

end