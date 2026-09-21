function [H_k] = get_measurement_model(r_k)
% Enum and constants
enum_dyn_model  = enum_dynamic_model();

switch r_k
    case enum_dyn_model.cv
        H_k = [1   0   0   0 ;
               0   0   1   0];

    case {enum_dyn_model.ct_r, enum_dyn_model.ct_l}
        H_k = [1   0   0   0   0 ;
               0   0   1   0   0];
end