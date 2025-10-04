function [Rho] = LDPC_opt_RHO(Objective,Inequality,Equality,Rh0,Lambda,InTer)
                              
% InTer=InTer+10;
Rho = Rh0;
beta = 0.9; % this parameter is used in the line search
Eps = 1e-1; % this parameter is used in the line search.
            % ...It is the sphere constraint. Can be higher
            % ...at the beginning and then get shorter.

% [~,D_v] = size(Lambda);
[~,D_c] = size(Rho);
% R_out = rate_calculation(Lambda,Rho);

n=0;
last_n = n;

% [~,Objective_grad,~] = feval(Objective, Rho);
[Equality_value] = feval(Equality,Rho);
[Inequality_value, ~] = feval(Inequality,Rho);

u = ones(D_c-1,1);  % this is the traspose of matrix A of the equality constraints, which is a vector here.
T = length(Inequality_value);
TT = length(Equality_value);
% iter = 0;
y = ones(T,1);   %initializing the vectors 'y' and 'z'
z = ones(TT,1);
gamma = -y'*Inequality_value;   % surrogate duality gap

%% reporting the parameters

% while (gamma >= tol || norm(Objective_grad+(Inequality_grad')*y+u*z)>=tol_feas || norm(1-u'*(Rho(2:end))') >= tol_feas)  % 3 stopping conditions are used
%     iter    
    
for interIter=1:InTer
    [~,Objective_grad,Objective_hess] = feval(Objective, Rho);
    [Equality_value] = feval(Equality,Rho);
    [Inequality_value, Inequality_grad] = feval(Inequality,Rho);
    
    R_out = rate_calculation(Lambda,Rho);   
    
    T = length(Inequality_value);
    
    if last_n > n
        temp_y = zeros(T,1);
        temp_y(1:D_c+n-1) = y(1:D_c+n-1);
        temp_y(D_c+n:end) = y(D_c+last_n:D_c+last_n+n-1);       
    else
        temp_y = ones(T,1);
        temp_y(1:D_c+last_n-1) = y(1:D_c+last_n-1);
        temp_y(D_c+n:D_c+n+last_n-1) = y(D_c+last_n:end);
    end
    y = temp_y;
    last_n = n;
    miu = 10 * T / gamma;
    
    coeff_mat = [Objective_hess          Inequality_grad'       u ;                         ...
                 diag(y)*Inequality_grad diag(Inequality_value) zeros(T,size(u,2));         ...
                 u'                      zeros(size(u',1),T)    zeros(size(u',1),size(u,2))];
             
    const_mat = [-Objective_grad-(Inequality_grad')*y-u*z;              ...
                 -diag(y)*Inequality_value-ones(T,1)/miu;               ...
                 1-u'*(Rho(2:end))'];
             
    direction = (coeff_mat)^(-1) * const_mat; %calculating the direction using the Newton's method
    delta_Rho = direction(1:D_c-1);
    delta_y = direction(D_c:D_c+T-1);
    delta_z = direction(D_c+T:end);
    s = Eps;                  %initializing the parameter for the line search   
    Rh0 = Rho';
    Rh0(2:end) = (Rho(2:end))' + s * delta_Rho;    % Lambda is a row vector, but delta_Lambda is a column vector
    Rh0 = (Rh0)';
    
    [Inequality_value, Inequality_grad] = feval(Inequality,Rho);
    
    y0 = y + s * delta_y;
    
    % doing line search to maintain feasibility and maintaining y>
    while( sum(Inequality_value(1:D_c-1)<=0) ~= D_c-1 || sum(y0>=0) ~= T)        
        s = s * beta ;     
        Rh0 = Rho';
        Rh0(2:end) = (Rho(2:end))' + s * delta_Rho;
        Rh0 = (Rh0)';        
        y0 = y + s * delta_y;
        [Inequality_value, Inequality_grad] = feval(Inequality,Rho);

    end
    
    Rho = Rh0;         %updating 'Lambda' and 'y' and 'z'
%     save(num2str(iter));
    y = y0;
    z = z + s * delta_z;
                          
    gamma = -y'*Inequality_value;
    
% 	n=0;
end
end