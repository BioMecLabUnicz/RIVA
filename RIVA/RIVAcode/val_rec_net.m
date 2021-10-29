function A_id_net_cell= val_rec_net(Normal_A_id,num_criteria)
% Given the normalised estimated matrix (Normal_A_id), the function 
% computes an array of matrices (A_id_net_cell), where, for each matirx, 
% new edges (i.e. entries) are added step-by-step from those with high 
% probability (by sorting the entries of the normalized matrix in 
% descending order by the absolute magnitude, the scalar input 
% num_criteria=1 for our aims): the matrix in the first cell corresponds to 
% the identified network with only one edge (the one with highest 
% probability), the matrix in the last cell corresponds to the full 
% identified network.

n_nodes=size(Normal_A_id,1);

for idx2=1:n_nodes
    Normal_A_id(idx2,idx2)=0;
end

if num_criteria==1
    
    
    %inc_coef2add=1;%n_nodes/2;
    if n_nodes<=30
        inc_coef2add=1;
        A_id_net_cell{(n_nodes-1)*(n_nodes/inc_coef2add)}=[];
    elseif n_nodes==100
        inc_coef2add=[1 2 2 3 3 4:35 40 50 75 100 200 500 800 1000 1500 2000 3000];
            
        A_id_net_cell{length(inc_coef2add)}=[];
    elseif n_nodes==1565
        inc_coef2add=[1:10, 15:5:50, 60:10:100 150:50:500 600:100:1000 2000:1000:10000 20000:10000:100000 200000:100000:500000 446345];

        A_id_net_cell{length(inc_coef2add)}=[];
    else
        
        inc_coef2add=input('vector containg the number of links to be added - strarting from the high-confidence predictions'); 
        A_id_net_cell{length(inc_coef2add)}=[];
    end
        
    % Compute the edge ranking matrix (Constr_list)
    [idx_row, idx_col, val_coef]=find(Normal_A_id);
    [val_coef_sort, idx_sort]=sort(val_coef,'descend');
    Constr_list=[idx_row(idx_sort) idx_col(idx_sort) val_coef_sort];
    num_coef2add=0;
    A_id_net=eye(n_nodes);
    for num_net_id=1:length(A_id_net_cell)
        if n_nodes<=30
            num_coef2add=num_coef2add+inc_coef2add;
        else       
            num_coef2add=num_coef2add+inc_coef2add(num_net_id);
        end
        if num_coef2add>size(Constr_list,1)
            num_coef2add=size(Constr_list,1);
        end
        for idx1=1:num_coef2add
            A_id_net(Constr_list(idx1,1),Constr_list(idx1,2))=Constr_list(idx1,3);
        end
        A_id_net_cell{num_net_id}=A_id_net;
    end
    
else
     A_id_net_cell{(n_nodes-1)}=[];
     Indx_mtx=zeros(n_nodes,n_nodes);
     for idx1=1:n_nodes
         [~, ind]=sort(Normal_A_id(idx1,:),'descend');
         Indx_mtx(idx1,:)=ind;
     end
     A_id_net=eye(n_nodes);
     for num_net_id=1:n_nodes-1
     
         idx_col=Indx_mtx(:,num_net_id);
         for idx1=1:n_nodes
             A_id_net(idx1,idx_col(idx1))=Normal_A_id(idx1,idx_col(idx1));
         end
         A_id_net_cell{num_net_id}=A_id_net;
     end
end
