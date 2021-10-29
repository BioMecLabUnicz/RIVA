
function [AUPR_m, AUROC_m, PPV, Se, Sp, FP]=val_aupr_auroc(original_network,A_id_cell,tag_norm)
% Given the original connectivity matrix (original_network matrix) and
% the reconstruced connectivity matrix (A_id_cell matrix), the function 
% evaluates the AUPR and AUROC by computing the PPV, Sensitivity, 
% Specificity and False Positives. The input tag_norm determines the way to
% normalizate the reconstructed connectivity matrix (= 1 for our aims, i.e. 
% normalization of each coefficient with respect to the values on
% the same row and column)

n_nodes=size(original_network,1);
TN=(n_nodes-1)*n_nodes-nnz(original_network);

% normalised estimated adjacency matrix 
Normal_A_id = norm_rec_mtx(A_id_cell,tag_norm);

% computes an array of matrices (A_id_net_cell), where, for each matrix, 
% new edges (i.e. entries) are added step-by-step from those with high 
% probability (by sorting the entries of the normalized matrix in 
% descending order by the absolute magnitude): the matrix in the first cell 
% corresponds to the identified network with only one edge (the one with 
% highest probability); the matrix in the last cell corresponds to the full 
% identified network.
A_id_net_cell= val_rec_net(Normal_A_id,1);

% performance indexes
PPV=zeros(1,size(A_id_net_cell,2));
Se=zeros(1,size(A_id_net_cell,2));
Sp=zeros(1,size(A_id_net_cell,2));
FP=zeros(1,size(A_id_net_cell,2));
for idx=1:size(A_id_net_cell,2)
    
    % compute reconstruction performance indexes:
    % PPV (positive predictive value), sensitivity, specificity and false
    % positives 
    if n_nodes==1565 % for big network compute the perfomance only for 
                     %  DIRECTED NETWORK
        Perf_idx= ppv_sens_spec_dir(original_network,A_id_net_cell{idx});

        PPV(idx)=Perf_idx(1);
        Se(idx)=Perf_idx(2);
        Sp(idx)=Perf_idx(3);
        FP(idx)=Perf_idx(4);
    else
        % compute the perfomance for UNDIRECTED, DIRECTED and SIGNED 
        % NETWORK, however for our aims we assume only DIRECTED NETWORK
        [Perf_idx, Res]= ppv_sens_spec(original_network,A_id_net_cell{idx});

        PPV(idx)=Perf_idx(2,1);
        Se(idx)=Perf_idx(2,2);
        Sp(idx)=Perf_idx(2,3);
        FP(idx)=size(Res.Dir.False_pos,1);
    end
end

% compute aupr and auroc indexes
AUPR_m=eval_aupr(PPV,Se);
AUROC_m=eval_auroc(Se,FP,TN);


