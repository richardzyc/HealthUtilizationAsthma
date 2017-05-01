function [b_new, pi_new] = Multinom_Est_Step(b, Z, Expl_mm_rev, pi, maxit, Overall_Iter)
K = size(pi,2);
param_values = b;
delta = ones(size(b));  %% Safety region paratemeter for each MN logistic regression coefficient.
mm_b = Expl_mm_rev*b.';
exp_mm_b = exp(mm_b);
log_lik = sum(sum(Z(:,1:K-1).*mm_b,2) - log(ones([size(Z,1),1]) + sum(exp_mm_b,2)), 1);
Z_pi_dif = Z(:,1:K-1) - pi(:,1:K-1);

differences = [];
dif = Inf;
MN_iter = 0;

while MN_iter < maxit && dif > 1
    MN_iter = MN_iter + 1;
    ['Overall Iteration: ',num2str(Overall_Iter),', Running MN Iteration: ',num2str(MN_iter)]
  
    %% Re-estimating each MN regression coefficient.
    for k = 1:K-1
        for j = 1:9
            v_kj = (sum(Z_pi_dif(:,k).*Expl_mm_rev(:,j)))/(sum(pi(:,k).*(ones(size(pi,1),1)-pi(:,k)).*Expl_mm_rev(:,j).*Expl_mm_rev(:,j)));
            abs_b_delta = min(abs(v_kj),delta(k,j));
            b_delta = abs_b_delta*sign(v_kj);
            delta(k,j) = max(2*abs_b_delta.',delta(k,j)/2);
            b(k,j) = b(k,j) + b_delta;
            mm_b = Expl_mm_rev*b.';
            exp_mm_b = exp(mm_b);
            pi_not_normed = [exp_mm_b, ones([size(exp_mm_b,1),1])];
            rowsums = sum(pi_not_normed,2);
            pi = zeros(size(pi_not_normed));
            for k_b = 1:K
                pi(:,k_b) = pi_not_normed(:,k_b)./rowsums;
            end
            Z_pi_dif = Z(:,1:K-1) - pi(:,1:K-1);
        end
    end
    cur_b = b;
    b
    delta
    
    cur_log_lik = sum(sum(Z(:,1:K-1).*mm_b,2) - log(ones([size(Z,1),1]) + sum(exp_mm_b,2)), 1);
    log_lik = [log_lik; cur_log_lik];
    dif = log_lik(size(log_lik, 1)) - log_lik(size(log_lik, 1) - 1);
    ['Last Multinomial Log Likelihood Improvement: ',num2str(dif)]
    differences = [differences; dif];
    param_values = [param_values; b];
end

b_new = b;
pi_new = pi;
end
