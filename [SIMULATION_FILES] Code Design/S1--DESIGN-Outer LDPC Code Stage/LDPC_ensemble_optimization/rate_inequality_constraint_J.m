function [ineq, ineq_grad]...
    = rate_inequality_constraint_J(lambda,                                ...
                                   MAIN_f_exit,                           ...
                                   MAIN_final_exit,MAIN_p,                ...
                                   WRTP_f_exit,                           ...
                                   WRTP_final_exit,WRTP_p)
                             
    D_v = length(lambda);
    n = length(MAIN_p)+length(WRTP_p);
    ineq = zeros(D_v+2*n-1,1); 
    ineq_grad = zeros(D_v+2*n-1,D_v-1);

    ineq(1:D_v-1) = -lambda(2:end);
    
    S = [MAIN_final_exit;WRTP_final_exit];
    p = [MAIN_p;WRTP_p];
    f_exit = [MAIN_f_exit;WRTP_f_exit];
    
    ineq(D_v:D_v+n-1) = S - p + 10^-10;
    ineq(D_v+n:end) = -S + p/exp(2) + 10^-10;  % convexity constraint
    ineq_grad((1:D_v-1),:) = -eye(D_v-1);
    ineq_grad((D_v:D_v+n-1),:) = f_exit(:,2:end);
    ineq_grad((D_v+n:end),:) = -f_exit(:,2:end);
end