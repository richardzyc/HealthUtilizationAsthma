function log_lik_overall = Overall_Log_Lik(pi_new, Cont_mm_rev, beta_new, surv_rev, Z_new, ic, Overall_Iter)
['Overall Iteration: ',num2str(Overall_Iter),', Calculating Overall Likelihood']

Z_new_for_PH = Z_new(ic,:);

R = max(surv_rev(:,10));
tau = surv_rev(:,3);
tau = tau/365;
tau = [tau, tau, tau, tau, tau, tau];
event_ind=surv_rev(:,4:9);

beta1 = beta_new(:,[1,6:20]);
exp_mm_beta1 = exp(Cont_mm_rev*beta1.');
mm_beta1 = log(exp_mm_beta1);
Z_1 = [Z_new_for_PH(:,1), Z_new_for_PH(:,1), Z_new_for_PH(:,1), Z_new_for_PH(:,1), Z_new_for_PH(:,1), Z_new_for_PH(:,1)];
log_lik1 = Z_1.*event_ind.*mm_beta1 - Z_1.*tau.*exp_mm_beta1;
log_Lik1 = sum(log_lik1,2);
beta2 = beta_new(:,[2,6:20]);
exp_mm_beta2 = exp(Cont_mm_rev*beta2.');
mm_beta2 = log(exp_mm_beta2);
Z_2 = [Z_new_for_PH(:,2), Z_new_for_PH(:,2), Z_new_for_PH(:,2), Z_new_for_PH(:,2), Z_new_for_PH(:,2), Z_new_for_PH(:,2)];
log_lik2 = Z_2.*event_ind.*mm_beta2 - Z_2.*tau.*exp_mm_beta2;
log_Lik2 = sum(log_lik2,2);
beta3 = beta_new(:,[3,6:20]);
exp_mm_beta3 = exp(Cont_mm_rev*beta3.');
mm_beta3 = log(exp_mm_beta3);
Z_3 = [Z_new_for_PH(:,3), Z_new_for_PH(:,3), Z_new_for_PH(:,3), Z_new_for_PH(:,3), Z_new_for_PH(:,3), Z_new_for_PH(:,3)];
log_lik3 = Z_3.*event_ind.*mm_beta3 - Z_3.*tau.*exp_mm_beta3;
log_Lik3 = sum(log_lik3,2);
beta4 = beta_new(:,[4,6:20]);
exp_mm_beta4 = exp(Cont_mm_rev*beta4.');
mm_beta4 = log(exp_mm_beta4);
Z_4 = [Z_new_for_PH(:,4), Z_new_for_PH(:,4), Z_new_for_PH(:,4), Z_new_for_PH(:,4), Z_new_for_PH(:,4), Z_new_for_PH(:,4)];
log_lik4 = Z_4.*event_ind.*mm_beta4 - Z_4.*tau.*exp_mm_beta4;
log_Lik4 = sum(log_lik4,2);
beta5 = beta_new(:,[5,6:20]);
exp_mm_beta5 = exp(Cont_mm_rev*beta5.');
mm_beta5 = log(exp_mm_beta5);
Z_5 = [Z_new_for_PH(:,5), Z_new_for_PH(:,5), Z_new_for_PH(:,5), Z_new_for_PH(:,5), Z_new_for_PH(:,5), Z_new_for_PH(:,5)];
log_lik5 = Z_5.*event_ind.*mm_beta5 - Z_5.*tau.*exp_mm_beta5;
log_Lik5 = sum(log_lik5,2);


log_Lik = [log_Lik1, log_Lik2, log_Lik3, log_Lik4, log_Lik5,  surv_rev(:,10)];
counts = hist(log_Lik(:,6),unique(log_Lik(:,6)));
log_Lik_array = mat2cell(log_Lik,counts,size(log_Lik,2));

surv_rev_log_lik = zeros(size(pi_new));

for r = 1:R
    surv_rev_log_lik(r,:) = sum(log_Lik_array{r}(:,1:5),1);
end

total_log_lik = Z_new.*log(pi_new)+surv_rev_log_lik;
log_lik_overall = sum(sum(total_log_lik,2));
end
