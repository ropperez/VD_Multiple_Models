function [Q] = get_process_noise(r_k,r_k_1,sigmaU,T,dim_nc,n_Q,ENUM)
% Enum
enum_dyn_model = ENUM.dynamic_model;

switch r_k
    case enum_dyn_model.cv
        switch r_k_1
            case enum_dyn_model.cv
                % Standard deviation
                sigma_cv = sigmaU.cv(n_Q);
                % Covariance matrix
                q        = [T^3/3   T^2/2 ;
                            T^2/2   T    ];
                I_dim    = eye(dim_nc);
                Q        = (sigma_cv)^2*kron(I_dim,q);

            case {enum_dyn_model.ct_r, enum_dyn_model.ct_l}
                % Standard deviation
                sigma_cv = sigmaU.cv(n_Q);
                % Covariance matrix
                q       = [T^3/3   T^2/2 ;
                           T^2/2   T    ];
                I_dim   = eye(dim_nc);
                Q       = (sigma_cv)^2*kron(I_dim,q);
        end

    case {enum_dyn_model.ct_r, enum_dyn_model.ct_l}
        % Standard deviation ct
        if r_k == enum_dyn_model.ct_r
            sigma_ct = sigmaU.ct_r(n_Q);
        elseif r_k == enum_dyn_model.ct_l
            sigma_ct = sigmaU.ct_l(n_Q);
        end

        switch r_k_1
            case enum_dyn_model.cv
                % Standard deviation
                sigma_cv = sigmaU.cv(n_Q);
                % Covariance matrix
                q        = [T^3/3   T^2/2 ;
                            T^2/2   T    ];
                I_dim    = eye(dim_nc);
                Q        = (sigma_cv)^2*blkdiag(kron(I_dim,q),(sigma_ct)^2*T);            

            case {enum_dyn_model.ct_r, enum_dyn_model.ct_l}
                % Standard deviation
                sigma_cv = sigmaU.cv(n_Q);
                % Covariance matrix
                q       = [T^3/3   T^2/2 ;
                           T^2/2   T    ];
                I_dim   = eye(dim_nc);
                Q       = (sigma_cv)^2*blkdiag(kron(I_dim,q),(sigma_ct)^2*T);
        end
end
end