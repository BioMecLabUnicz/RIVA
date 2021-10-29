function AUROC=eval_auroc(Se,FP,TN)
% compute the auroc index, the area under ROC curve, i.e. true positive 
% rate (Sensitivity) vs. false positive rate.

AUROC=0;

for idx=1:length(Se)
    if idx==1
       A_idx=Se(idx)*FP(idx)/TN;
    elseif Se(idx)==1 && FP(idx)==0
        AUROC=1;
        break
    else
        A_idx=Se(idx)*((FP(idx)-FP(idx-1))/TN);
    end
    AUROC=AUROC+A_idx;
end

