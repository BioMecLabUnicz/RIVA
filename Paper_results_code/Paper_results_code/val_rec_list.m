function A_id_ones=val_rec_list(A_rank,n_edges)

% create a matrix (A_id_ones) with entries equal to 1  
% corresponding to the n edges (n=n_edges) with high probability of being 
% true positives
[idx_row, idx_col, val_coef]=find(A_rank);
   
[val_coef_sort, idx_sort]=sort(val_coef,'descend');

n_nodes=size(A_rank,1);
A_id_ones=eye(n_nodes);

if n_edges>length(val_coef_sort) 
            
    n_edges=length(val_coef_sort);
end        
for idx=1:n_edges
    A_id_ones(idx_row(idx_sort(idx)),idx_col(idx_sort(idx)))=1;

end