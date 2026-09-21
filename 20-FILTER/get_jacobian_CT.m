function [J] = get_jacobian_CT(x_k,t,dim,dim_c)
% Compute coordinated turn jacobian for transition matrix
%   x_k = [x,x_dot,y,y_dot,w]
%   x_t  = x+vx*sin(w*T)/w-vy*(1-cos(w*T))/w;
%   vx_t = vy*cos(w*T)-vy*sin(w*T);
%   y_t  = vx*(1-cos(w*T))/w + y + vy*sin(w*T)/w ;
%   vy_t = vx*sin(w*T) +vx*cos(w*T);
%   r_t  = sqrt(x_t^2+y_t^2);
%   w_t  = sqrt(vx_t^2+vy_t^2)/r_t;

% Dimension
L_x_k = length(x_k);

% Extracción de variables del Vector de Estado en Frame R1P
x       = x_k(1);
x_dot   = x_k(2);
y       = x_k(1+dim_c);
y_dot   = x_k(2+dim_c);
w       = x_k(end);

% Jacobian initialisation
J       = zeros(L_x_k);

% Avoiding division by practically zero, w --> 0
if abs(w) > sqrt(eps)
    % x(t)
    J(1,1)       =                                                                             1; % d(x(t))/dx           
    J(1,2)       =                                                                    sin(t*w)/w; % d(x(t))/dx_dot         
    J(1,1+dim_c) =                                                                             0; % d(x(t))/dy         
    J(1,2+dim_c) =                                                              (cos(t*w) - 1)/w; % d(x(t))/dy_dot
    J(1,dim)     = ((w*t*cos(w*t) - sin(w*t))*x_dot + (1 - cos(w*t) - w*t*sin(w*t))*y_dot) / w^2; % d(x(t))/dw
    
    % x_dot(t)
    J(2,1)       =                                      0; % d(x_dot(t))/dx          
    J(2,2)       =                               cos(t*w); % d(x_dot(t))/dx_t         
    J(2,1+dim_c) =                                      0; % d(x_dot(t))/dy         
    J(2,2+dim_c) =                              -sin(t*w); % d(x_dot(t))/dy_dot
    J(2,dim)     = (-x_dot*sin(w*t) - y_dot*cos(w*t)) * t; % d(x_dot(t))/dw
    
    % y(t)
    J(1+dim_c,1)       =                                                                                     0; % d(y(t))/dx                
    J(1+dim_c,2)       =                                                                     -(cos(t*w) - 1)/w; % d(y(t))/dx_dot                        
    J(1+dim_c,1+dim_c) =                                                                                     1; % d(y(t))/dy               
    J(1+dim_c,2+dim_c) =                                                                            sin(t*w)/w; % d(y(t))/dy_dot                  
    J(1+dim_c,dim)     = (w*t*(x_dot*sin(w*t) + y_dot*cos(w*t)) - (x_dot*(1-cos(w*t)) + y_dot*sin(w*t))) / w^2; % d(y(t))/dw                   
    
    % y_dot(t)
    J(2+dim_c,1)       =                                     0; % d(y_dot)/dx            
    J(2+dim_c,2)       =                              sin(t*w); % d(y_dot)/dx_dot           
    J(2+dim_c,1+dim_c) =                                     0; % d(y_dot)/dy           
    J(2+dim_c,2+dim_c) =                              cos(t*w); % d(y_dot)/dy_dot        
    J(2+dim_c,dim)     = (x_dot*cos(w*t) - y_dot*sin(w*t)) * t; % d(y_dot)/dw       
    
    % w(t)
    J(5,1)       = 0; % d(w(t))/dx            
    J(5,2)       = 0; % d(w(t))/dx_dot           
    J(5,1+dim_c) = 0; % d(w(t))/dy  
    J(5,2+dim_c) = 0; % d(w(t))/dy_dot        
    J(5,dim)     = 1; % d(w(t))/dw       

else
    % Avoiding division by practically zero, take the limit of w --> 0
    % L'Hopital's rule for the limit

    % x(t)
    J(1,1)       =              1; % d(x(t))/dx           
    J(1,2)       =              t; % d(x(t))/dx_dot         
    J(1,1+dim_c) =              0; % d(x(t))/dy         
    J(1,2+dim_c) =              0; % d(x(t))/dy_dot
    J(1,dim)     = -0.5 * y_dot*t^2; % d(x(t))/dw
    
    % x_dot(t)
    J(2,1)       =          0; % d(x_dot(t))/dx          
    J(2,2)       =          1; % d(x_dot(t))/dx_t         
    J(2,1+dim_c) =          0; % d(x_dot(t))/dy         
    J(2,2+dim_c) =          0; % d(x_dot(t))/dy_dot
    J(2,dim)     = -y_dot * t; % d(x_dot(t))/dw
    
    % y(t)
    J(1+dim_c,1)       =               0; % d(y(t))/dx                
    J(1+dim_c,2)       =               0; % d(y(t))/dx_dot                        
    J(1+dim_c,1+dim_c) =               1; % d(y(t))/dy               
    J(1+dim_c,2+dim_c) =               t; % d(y(t))/dy_dot                  
    J(1+dim_c,dim)     = 0.5 * x_dot*t^2; % d(y(t))/dw                   
    
    % y_dot(t)
    J(2+dim_c,1)       =         0; % d(y_dot)/dx            
    J(2+dim_c,2)       =         0; % d(y_dot)/dx_dot           
    J(2+dim_c,1+dim_c) =         0; % d(y_dot)/dy           
    J(2+dim_c,2+dim_c) =         1; % d(y_dot)/dy_dot        
    J(2+dim_c,dim)     = x_dot * t; % d(y_dot)/dw       
    
    % w(t)
    J(dim,1)       = 0; % d(w(t))/dx            
    J(dim,2)       = 0; % d(w(t))/dx_dot           
    J(dim,1+dim_c) = 0; % d(w(t))/dy  
    J(dim,2+dim_c) = 0; % d(w(t))/dy_dot        
    J(dim,5)       = 1; % d(w(t))/dw  
end

end