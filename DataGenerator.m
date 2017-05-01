n=3000;% number of people

CONT=csvread('Control Variables.csv',1,0);
EXPL=csvread('Explanatory Variables.csv',1,0);
SURV=csvread('Survival Observations.csv',1,0);
sa=mnrnd(n, CONT(:,11)/sum(CONT(:,11)));
loc=EXPL(:,9)/sum(EXPL(:,9));
eventProb=[];
[r,c]=size(SURV);
for i=1:r
    eventProb=[eventProb;[SURV(i,1:7)/sum(SURV(i,1:7))]];
end
surv=[];
Cont_mm=[];
Expl_mm_rev=[];
pat=0;
for i=1:length(sa)
    
    for j=1:sa(i)
        pat=pat+1;
        Expl_mm_rev=[Expl_mm_rev; [EXPL(find(mnrnd(1,loc)),1:8), normrnd(EXPL(find(mnrnd(1,loc)),10), EXPL(find(mnrnd(1,loc)),11))]];
        prev=[0,0,0,0,0,0];
        year=2015;
        day=0;
        for z=1:round(SURV(i,10))
            daysBetween=exprnd(SURV(i,9));
            if day+daysBetween>365
                year=year+1;
                day=day+daysBetween-365;
            else
                day=day+daysBetween;
            end
            Cont_mm=[Cont_mm;[CONT(i,1:8),prev,CONT(i,9:10)]];
            prev=mnrnd(1,eventProb(i,:));
            prev=prev(1:6);
            surv=[surv;[[pat, year, daysBetween, prev,pat]]];
        end
    end
    
end
Expl_mm_rev(Expl_mm_rev(:,9)<0,9)=0.01;
Expl_mm_rev(Expl_mm_rev(:,9)>1,9)=1;