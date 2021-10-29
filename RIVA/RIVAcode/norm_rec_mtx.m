function Normal_A_id = norm_rec_mtx(A_id,tag_norm)
% Given the matrix A_id, the function returns the normalised estimated 
% adjacency matrix named Normal_A_id. The input tag_norm is a scalar that 
% determines the way to normalize the matrix A_id: 1, i.e. normalization of 
% each coefficient with respect to the values on the same row and column; 2
% normalization of each coefficient with respect to the values on the same 
% row; 3, i.e. the abs of each coefficient. For our aims, tag_norm=1.

n_nodes=size(A_id,1);
if tag_norm==1
    
    % Compute the normalized A matrix
    Normal_A_id=zeros(n_nodes,n_nodes);
    for idx1=1:n_nodes
        row_norm=norm(A_id(idx1,:));
        for idx2=1:n_nodes
            % Normalization of each coefficient with respect to the values on
            % the same row and column
            Normal_A_id(idx1,idx2)=abs(A_id(idx1,idx2))/sqrt(row_norm*norm(A_id(:,idx2)));
        end
    end

%     Normal_A_id=zeros(n_nodes,n_nodes);
%     for idx1=1:n_nodes
%         % Normalization of each coefficient with respect to the values on
%         % the same row
%         Normal_A_id(idx1,:)=abs(A_id(idx1,:))/norm(A_id(idx1,:));
%     end
%     for idx1=1:n_nodes
%         % Normalization of each coefficient with respect to the values on
%         % the same column
%         Normal_A_id(:,idx1)=Normal_A_id(:,idx1)/norm(A_id(:,idx1));
%     end
elseif  tag_norm==2
    
    % Compute the normalized A matrix
    Normal_A_id=zeros(n_nodes,n_nodes);
    for idx1=1:n_nodes
        % Normalization of each coefficient with respect to the values on
        % the same row
        Normal_A_id(idx1,:)=abs(A_id(idx1,:))/norm(A_id(idx1,:));
    end
else
    Normal_A_id=abs(A_id);
end


         

