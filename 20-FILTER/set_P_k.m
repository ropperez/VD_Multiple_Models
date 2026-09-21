function [P_k_dim] = set_P_k(P_k,r_k,dim)
enum_dyn_model = enum_dynamic_model();

P_k_dim = zeros(dim);

if dim == 7 || dim == 6
    if r_k == enum_dyn_model.ca
        P_k_dim(1:6,1:6) = P_k;
    elseif r_k == enum_dyn_model.cv || r_k == enum_dyn_model.ct_l || r_k == enum_dyn_model.ct_r
        P_k_dim(1:2,1:2) = P_k(1:2,1:2);
        P_k_dim(4:5,4:5) = P_k(3:4,3:4);
    end
elseif dim == 5
    if r_k == enum_dyn_model.cv
        P_k_dim(1:4,1:4) = P_k(1:4,1:4);
    elseif r_k == enum_dyn_model.ct_l || r_k == enum_dyn_model.ct_r
        P_k_dim = P_k;
    end    
end

end
