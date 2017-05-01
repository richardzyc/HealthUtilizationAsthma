function [param_values, differences, log_lik] = Initial_PH_Est_Step(surv, Cont_mm, maxit)
tau = surv(:,3); %% Column of interarrival times (in days)
tau = tau/365; %% Changing unit of interarrival times from days to years
event_ind=surv(:,4:9); %% Columns of event indicators of 1's and 0's.

param_values = [zeros([6, 16])]; %% 6 x 16 matrix of 0's containing PH Coefficients, one row for each event type and 16 control parameters.
beta = param_values;
delta = ones([6, 16]); %% The safety region parameter for each PH parameter as described in Algorithm 1.
mm_beta = Cont_mm*beta.';
exp_mm_beta = exp(mm_beta);
log_lik =  sum(event_ind.*mm_beta -([tau, tau, tau, tau, tau, tau].*exp_mm_beta));
differences = [];
dif = Inf;
PH_iter = 0;
cur_beta = beta;
cur_delta = delta;
old_log_lik = log_lik;

while PH_iter < maxit && (sum(dif) > .25 || max(max(delta)) > .01)
    PH_iter = PH_iter + 1;
    ['Running Initial PH Iteration: ',num2str(PH_iter)]
    % Calculating Parameter Updates
    for s = 1:6
        mm_beta_s = Cont_mm * beta(s,:).';
        exp_mm_beta_s = exp(mm_beta_s);
        for j = 1:16
            Cont_mm_j = Cont_mm(:,j);
            v_js = - ((sum(event_ind(:,s).*Cont_mm_j) - sum(tau.*Cont_mm_j.*exp_mm_beta_s))./(-sum(tau.*(Cont_mm_j.^2).*exp_mm_beta_s)));
            abs_beta_delta = min(abs(v_js),delta(s,j));
            beta_delta = abs_beta_delta * sign(v_js);
            delta(s,j) = max(2*abs_beta_delta,delta(s,j)/2);
            beta(s,j) = beta(s,j) + beta_delta;
            mm_beta_s = Cont_mm * beta(s,:).';
            exp_mm_beta_s = exp(mm_beta_s);
            cur_log_lik(:,s) = sum(event_ind(:,s).*mm_beta_s -(tau.*exp_mm_beta_s));
            if (cur_log_lik(:,s) - old_log_lik(:,s)) < 0
                beta(s,j) = cur_beta(s,j);
                delta(s,j) = cur_delta(s,j)/2;
                cur_log_lik(:,s) = old_log_lik(:,s);
                mm_beta_s = Cont_mm * beta(s,:).';
                exp_mm_beta_s = exp(mm_beta_s);
            end
            old_log_lik(:,s) = cur_log_lik(:,s);
        end
    end
    cur_beta = beta;
    cur_delta = delta;
    beta
    delta

    % Re-calculating Log Likelihood & Resetting Params
    log_lik = [log_lik; cur_log_lik];
    dif = log_lik(size(log_lik, 1),:) - log_lik(size(log_lik, 1) - 1,:);
    old_log_lik = log_lik(size(log_lik, 1),:);
    ['Last Log Likelihood Improvement: ',num2str(dif)]
    differences = [differences; dif];
    param_values = [param_values; cur_beta];
end

end

