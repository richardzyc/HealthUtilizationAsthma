%% Overall Script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% surv contains the [ Patient_ID, Year, Interarrival_Times, CL_Indicator, ER_Indicator, HO_Indicator, RX_Indicator, NP_Indicator, PO_Indicator, St_Patient_ID]
%% Cont_mm is a (Number of Lifetimes x 16) matrix containing the [Baseline_Indicator, Age_Grp_1_Indicator, Age_Grp_2_Indicator, Black_Indicator, Other_Race_Indicator, Healthy_Indicator, Minor_Ill_Indicator, Severe_Ill_Indicator, Blind_Disabled_Indicator, Foster_Care_Indicator, Prev_CL_Event, Prev_ER_Event, Prev_HO_Event, Prev_RX_Event, Prev_NP_Event, Prev_PO_Event]
%% Expl_mm_rev is a (R x 9) matrix containing the Baseline_Indicator, LA_Indicator, MS_Indicator, MN_Indicator, NC_Indicator, TN_Indicator, Suburban_Indicator, Rural_Indicator, Travel_Distance]

%% Run Initial PH Parameter Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[param_values, differences, log_lik] = Initial_PH_Est_Step(surv, Cont_mm, 100)
param_inits = param_values((((size(param_values,1)/6)-1)*6 + 1):size(param_values,1),:);

for K = 5
    for N = 1:2
        %% Initializing Cluster Membership Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        id_st = surv(:,10);  %% The 10th column of of surv is a unique indicator for each patient/state combination. There are R of these.
        [id_st_un, ia, ic] = unique(id_st, 'rows', 'stable');
        R = max(id_st_un);

	%% Here, we randomly assign each patient to a cluster.
        Z = zeros([size(ia,1),K]);
        for r = 1:R
            Z(r,:) = randperm(K); 
        end
        Z(Z > 1) = 0;
        save(['Runs_K_',num2str(K),'\N_',num2str(N),'\Z_init.mat'],'Z');
        Z_for_PH = Z(ic,:);

	%% We set the initial baseline PH parameter the same across all clusters here.
        Cont_mm_for_PH = [repmat(Cont_mm(:,1),1,K), Cont_mm(:,2:16)]; 
        beta = [repmat(param_inits(:,1),1,K), param_inits(:,2:16)]
        delta = ones([6, (K+15)]) %% The safety region parameter for each PH parameter.

        Overall_Iter = 0;
	
	%% Here, we set the multinomial logistic regression coefficients to 0.
        b = zeros([(K-1), 9]);
        pi = 1/K * ones(size(Z));

        log_lik_overall =  Overall_Log_Lik(pi, Cont_mm, beta, surv, Z, ic, Overall_Iter)

        %% EM Algorithm: 1st Iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Overall_Iter = Overall_Iter + 1;

        beta_new = PH_Est_Step_K_5(surv, beta, delta, Cont_mm_for_PH, Z_for_PH, 5, Overall_Iter);
        save(['Runs_K_',num2str(K),'\N_',num2str(N),'\beta_',num2str(Overall_Iter),'.mat'],'beta_new');

        [b_new, pi_new] = Multinom_Est_Step(b,Z, Expl_mm_rev, pi, 100, Overall_Iter);
        save(['Runs_K_',num2str(K),'\N_',num2str(N),'\b_',num2str(Overall_Iter),'.mat'],'b_new');
        save(['Runs_K_',num2str(K),'\N_',num2str(N),'\pi_',num2str(Overall_Iter),'.mat'],'pi_new');

        Z_new = Estep(Cont_mm, surv, beta_new, pi_new, Overall_Iter);
        save(['Runs_K_',num2str(K),'\N_',num2str(N),'\Z_',num2str(Overall_Iter),'.mat'],'Z_new');
        Z_new_for_PH = Z_new(ic,:);

        cur_log_lik_overall = Overall_Log_Lik(pi_new, Cont_mm, beta_new, surv, Z_new, ic, Overall_Iter);
        log_lik_overall = [log_lik_overall; cur_log_lik_overall];

        difference = cur_log_lik_overall - log_lik_overall(size(log_lik_overall,1)-1);
        ['Likelihood Improvement = ',num2str(difference/abs(cur_log_lik_overall)*100),'%']

        %% EM Algorithm: Subsequent Iterations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while Overall_Iter < 100 && (Overall_Iter < 5 || difference/abs(cur_log_lik_overall) > 1e-10) %% thresh determines how tight the convergence criterion should be of the overall (complete) likelihood function. That is, once the increase in likelihood is less than "thresh", the algorithm stops.

            Overall_Iter = Overall_Iter + 1;

            beta_new = PH_Est_Step_K_5(surv, beta_new, delta, Cont_mm_for_PH, Z_new_for_PH, 5, Overall_Iter);
            save(['Runs_K_',num2str(K),'\N_',num2str(N),'\beta_',num2str(Overall_Iter),'.mat'],'beta_new');

            [b_new, pi_new] = Multinom_Est_Step(b_new, Z_new, Expl_mm_rev, pi_new, 100, Overall_Iter);
            save(['Runs_K_',num2str(K),'\N_',num2str(N),'\b_',num2str(Overall_Iter),'.mat'],'b_new');
            save(['Runs_K_',num2str(K),'\N_',num2str(N),'\pi_',num2str(Overall_Iter),'.mat'],'pi_new');

            Z_new = Estep(Cont_mm, surv, beta_new, pi_new, Overall_Iter);
            save(['Runs_K_',num2str(K),'\N_',num2str(N),'\Z_',num2str(Overall_Iter),'.mat'],'Z_new');
            Z_new_for_PH = Z_new(ic,:);

            cur_log_lik_overall = Overall_Log_Lik(pi_new, Cont_mm, beta_new, surv, Z_new, ic, Overall_Iter);
            log_lik_overall = [log_lik_overall; cur_log_lik_overall];
            difference = cur_log_lik_overall - log_lik_overall(size(log_lik_overall,1)-1);
            ['Likelihood Improvement = ',num2str(difference/abs(cur_log_lik_overall)*100),'%']
        end

        save(['Runs_K_',num2str(K),'\N_',num2str(N),'\log_lik_overall.mat'], 'log_lik_overall');
    end
end