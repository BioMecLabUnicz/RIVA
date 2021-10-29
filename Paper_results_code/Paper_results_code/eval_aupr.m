function AUPR=eval_aupr(PPV,Se)
% compute the aupr index, the area under PPV vs Se

AUPR=0;

for idx=1:length(PPV)
    if idx==1
       A_idx=PPV(idx)*Se(idx);
    else
        A_idx=PPV(idx)*(Se(idx)-Se(idx-1));
    end
    AUPR=AUPR+A_idx;
end

