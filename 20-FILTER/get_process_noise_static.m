function [Q] = get_process_noise_static(r_k,sigmaU,T,dim_nc,n_Q,ENUM)
enum_dyn_model  = ENUM.dynamic_model;

switch r_k
    case enum_dyn_model.cv
        % Nearly constant velocity model:
        sigma_cv    = sigmaU.cv(n_Q);

        % Covariance matrix
        I_dim       = eye(dim_nc);

        q = [T^3/3   T^2/2 ;
             T^2/2   T    ];

        Q = (sigma_cv)^2*blkdiag(kron(I_dim,q),0);

    case {enum_dyn_model.ct_r,enum_dyn_model.ct_l}
        % Standard deviation
        sigma_cv = sigmaU.cv(n_Q);
        if r_k == enum_dyn_model.ct_r
            sigma_ct    = sigmaU.ct_r(n_Q);
        elseif r_k == enum_dyn_model.ct_l
            sigma_ct    = sigmaU.ct_l(n_Q);
        end

        % Covariance matrix
        I_dim       = eye(dim_nc);

        q = [T^3/3   T^2/2 ;
             T^2/2   T    ];

        Q = (sigma_cv)^2*blkdiag(kron(I_dim,q),(sigma_ct)^2*T);
end
end