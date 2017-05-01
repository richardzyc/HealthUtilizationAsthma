function [Z_new] = Estep(Cont_mm_rev, surv_rev, beta, pi, Overall_Iter)
['Overall Iteration: ',num2str(Overall_Iter),', Calculating new Z matrix']
K = size(pi,2);
surv_Lik = zeros(size(pi));
R = max(surv_rev(:,10));
tau = surv_rev(:,3);
tau = tau/365;
tau = [tau, tau, tau, tau, tau, tau];
event_ind=surv_rev(:,4:9);

Lik = [];
for k = 1:K
    betak = beta(:,[k,(K+1):(K+15)]);
    exp_mm_betak = exp(Cont_mm_rev*betak.');
    mm_betak = log(exp_mm_betak);
    lik_k = exp(event_ind.*mm_betak - tau.*exp_mm_betak);
    Lik_k = prod(lik_k,2);
    Lik = [Lik, Lik_k];
end

Lik = [Lik, surv_rev(:,10)];
counts = hist(Lik(:,K+1),unique(Lik(:,K+1)));
Likarray = mat2cell(Lik,counts,size(Lik,2));

for r = 1:R
    surv_Lik(r,:) = prod(Likarray{r}(:,1:K),1);
end

Z_raw = surv_Lik.*pi;
[row, col] = find(Z_raw == Inf);
Z_raw_v2 = Z_raw;
Z_raw_v2(row,:) = 0;
Z_raw_v2(find(Z_raw == Inf)) = 1;
rowsums = sum(Z_raw_v2,2);
Z_new = zeros(size(Z_raw_v2));
for k = 1:K
    Z_new(:,k) = Z_raw_v2(:,k)./rowsums;
end

end