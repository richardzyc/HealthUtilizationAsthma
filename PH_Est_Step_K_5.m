function beta_new = PH_Est_Step_K_5(surv_rev, beta, delta, Cont_mm_for_PH, Z_for_PH, maxit, Overall_Iter)

K = size(Z_for_PH,2);
tau = surv_rev(:,3);
tau = tau/365;
tau = [tau, tau, tau, tau, tau, tau];  %% Tau must be replicated for each event type (6 times).
event_ind=surv_rev(:,4:9);

log_lik = 0;
for k = 1:5
    mm_beta_k = Cont_mm_for_PH(:,[k, (K+1):(15+K)])*beta(:,[k, (K+1):(15+K)]).';
    exp_mm_beta_k = exp(mm_beta_k);
    Z_k = [Z_for_PH(:,k), Z_for_PH(:,k), Z_for_PH(:,k), Z_for_PH(:,k), Z_for_PH(:,k), Z_for_PH(:,k)];
    log_lik_k = sum(Z_k.*event_ind.*mm_beta_k -(tau.*Z_k.*exp_mm_beta_k));
    log_lik = log_lik + log_lik_k;
end
dif = Inf;
PH_iter = 0;

%% Like tau, we must replicate the estimate for Z^m 6 times, once for each event.
Z_1 = [Z_for_PH(:,1), Z_for_PH(:,1), Z_for_PH(:,1), Z_for_PH(:,1), Z_for_PH(:,1), Z_for_PH(:,1)];  
Z_2 = [Z_for_PH(:,2), Z_for_PH(:,2), Z_for_PH(:,2), Z_for_PH(:,2), Z_for_PH(:,2), Z_for_PH(:,2)];
Z_3 = [Z_for_PH(:,3), Z_for_PH(:,3), Z_for_PH(:,3), Z_for_PH(:,3), Z_for_PH(:,3), Z_for_PH(:,3)];
Z_4 = [Z_for_PH(:,4), Z_for_PH(:,4), Z_for_PH(:,4), Z_for_PH(:,4), Z_for_PH(:,4), Z_for_PH(:,4)];
Z_5 = [Z_for_PH(:,5), Z_for_PH(:,5), Z_for_PH(:,5), Z_for_PH(:,5), Z_for_PH(:,5), Z_for_PH(:,5)];

while PH_iter < (maxit) && sum(dif) > 1
    PH_iter = PH_iter + 1;
    ['Running Overall Iteration: ',num2str(Overall_Iter),', PH Iteration: ',num2str(PH_iter)]
    mm_beta_nk = Cont_mm_for_PH(:,[(K+1):(15+K)])*beta(:,[(K+1):(15+K)]).';
    
    %% Re-estimating baseline parameters.
    for k = 1:K
        mm_beta_0k = Cont_mm_for_PH(:,k)*beta(:,k).';
        mm_beta_k = mm_beta_nk + mm_beta_0k;
        exp_mm_beta_k = exp(mm_beta_k);
        Z_k = [Z_for_PH(:,k), Z_for_PH(:,k), Z_for_PH(:,k), Z_for_PH(:,k), Z_for_PH(:,k), Z_for_PH(:,k)];
        v_s = - ((sum(Z_k.*event_ind) - sum(Z_k.*tau.*exp_mm_beta_k))./(-sum(Z_k.*tau.*exp_mm_beta_k)));
        abs_beta_delta = min(abs(v_s).',delta(:,1));
        beta_delta = abs_beta_delta .* (sign(v_s).');
        beta(:,k) = beta(:,k) + beta_delta;
        delta(:,k) = max(2*abs_beta_delta, delta(:,k)/2);
    end
    
    mm_beta_01 = Cont_mm_for_PH(:,1)*beta(:,1).';
    exp_mm_beta_01 = exp(mm_beta_01);
    mm_beta_02 = Cont_mm_for_PH(:,2)*beta(:,2).';
    exp_mm_beta_02 = exp(mm_beta_02);
    mm_beta_03 = Cont_mm_for_PH(:,3)*beta(:,3).';
    exp_mm_beta_03 = exp(mm_beta_03);
    mm_beta_04 = Cont_mm_for_PH(:,4)*beta(:,4).';
    exp_mm_beta_04 = exp(mm_beta_04);
    mm_beta_05 = Cont_mm_for_PH(:,5)*beta(:,5).';
    exp_mm_beta_05 = exp(mm_beta_05);

    %% Re-estimating control coefficient PH parameters.
    for j = 1:15
        mm_beta_nk = Cont_mm_for_PH(:,[(K+1):(15+K)])*beta(:,[(K+1):(15+K)]).';
        exp_mm_beta_nk = exp(mm_beta_nk);
        exp_mm_beta_1 = exp_mm_beta_01.*exp_mm_beta_nk;
        exp_mm_beta_2 = exp_mm_beta_02.*exp_mm_beta_nk;
        exp_mm_beta_3 = exp_mm_beta_03.*exp_mm_beta_nk;
        exp_mm_beta_4 = exp_mm_beta_04.*exp_mm_beta_nk;
        exp_mm_beta_5 = exp_mm_beta_05.*exp_mm_beta_nk;
        exp_mm_beta = Z_1.*exp_mm_beta_1 + Z_2.*exp_mm_beta_2 + Z_3.*exp_mm_beta_3 + Z_4.*exp_mm_beta_4 + Z_5.*exp_mm_beta_5;
        Cont_mm_j = [Cont_mm_for_PH(:,j+K), Cont_mm_for_PH(:,j+K), Cont_mm_for_PH(:,j+K), Cont_mm_for_PH(:,j+K), Cont_mm_for_PH(:,j+K), Cont_mm_for_PH(:,j+K)];
        v_s = - ((sum(event_ind.*Cont_mm_j) - sum(tau.*Cont_mm_j.*exp_mm_beta))./(-sum(tau.*Cont_mm_j.*Cont_mm_j.*exp_mm_beta)));
        abs_beta_delta = min(abs(v_s).', delta(:, j+K));
        beta_delta = abs_beta_delta.*(sign(v_s).');
        beta(:,j+K) = beta(:, j+K) + beta_delta;
        delta(:,j+K) = max(2*abs_beta_delta, delta(:,j+K)/2);
    end
    beta
    delta

    % Re-calculating Log Likelihood
    new_log_lik = 0;
    for k = 1:K
        mm_beta_k = Cont_mm_for_PH(:,[k, (K+1):(15+K)])*beta(:,[k, (K+1):(15+K)]).';
        exp_mm_beta_k = exp(mm_beta_k);
        Z_k = [Z_for_PH(:,k), Z_for_PH(:,k), Z_for_PH(:,k), Z_for_PH(:,k), Z_for_PH(:,k), Z_for_PH(:,k)];
        log_lik_k = sum(Z_k.*event_ind.*mm_beta_k -(tau.*Z_k.*exp_mm_beta_k));
        new_log_lik = new_log_lik + log_lik_k;
    end
    log_lik = [log_lik; new_log_lik];
    dif = log_lik(size(log_lik, 1),:) - log_lik(size(log_lik, 1) - 1, :);
    ['Last PH Log Likelihood Improvement: ',num2str(dif)]
end
beta_new = beta;

end