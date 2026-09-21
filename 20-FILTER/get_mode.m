function [r_k] = get_mode(r_k_1,mu)
% Probability of transition
p_transition    = rand; 

% Motion model selection based on mode transition matrix
r_k = find(p_transition<= cumsum(mu(r_k_1,:)),1);
end