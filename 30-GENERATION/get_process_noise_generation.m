function [Q] = get_process_noise_generation(dyn_model,sigmaU,T,dim_nc,sim_r,N_q,ENUM)
enum_dyn_model  = ENUM.dynamic_model;

switch dyn_model
    case enum_dyn_model.cv
        % Nearly constant velocity model:
        sigma_cv    = sigmaU.cv(N_q);

        % Covariance matrix
        I_dim       = eye(dim_nc);

        q = [T^3/3   T^2/2 ;
             T^2/2   T    ];

        Q = (sigma_cv)^2*blkdiag(kron(I_dim,q),0);

    case {enum_dyn_model.ct_r,enum_dyn_model.ct_l}
        % Standard deviation
        sigma_cv = sigmaU.cv(N_q);
        if dyn_model == enum_dyn_model.ct_r
            sigma_ct    = sigmaU.ct_r(N_q);
        elseif dyn_model == enum_dyn_model.ct_l
            sigma_ct    = sigmaU.ct_l(N_q);
        end

        % Covariance matrix
        I_dim       = eye(dim_nc);

        q = [T^3/3   T^2/2 ;
             T^2/2   T    ];

        Q = (sigma_cv)^2*blkdiag(kron(I_dim,q),(sigma_ct)^2*T);
end
end
