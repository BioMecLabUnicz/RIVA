function [A_id,B_id,norm_est_err_dy,sol_found,Meas_mtx_sim,Err_mtx]...
    = riva(Meas_mtx,Meas_mtx_sim_pr,null_constr_list,num_tps_4exp)

% The function computes the connectivity matrix A and exogenous
% perturbation matrix B by using the RLS and IV techniques, with imposing  
% the desired zero entries of the connectivity matrix A and of the input 
% matrix  B contained in the list named null_constr_list (if no information  
% is available for A and B, then the list is empty).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT VARIABLES
%
% Meas_mtx:         Measurements matrix, each column contains the 
%                   measurements for a node
%
% Meas_mtx_sim_pr   Simulated measurements matrix, each column contains the 
%                   measurements for a node, obtained by the model
%                   identified previously, at first it equals to Meas_mtx
%                                 
% null_constr_list: List that contains the coefficients of the connectivity
%                   matrix A and of the input matrix B to be nullified  
%                   (emppty, if no information available)        
%  
% num_tps_4exp:     a vector containing the number of time points for each 
%                   experiment (time-series),  num_tps_4exp(i)= number of
%                   time points for the i-th time-series. The number of 
%                   rows of Meas_mtx is equal to sum(num_tps_4exp).
%
%
% OUTPUT VARIABLES
%
% A_id:            The identified connectivity matrix
%               
%
% B_id:            The identified input matrix   
%               
% norm_est_err_dy: Estimation error generated by the
%                  identified model 
%
% sol_found:       1 - the matrix A_id is stable, then compute the simulated
%                      measurements matrix (Meas_mtx_sim) by simulating the  
%                      time evolution of the identified system
%                  0 - the matrix A_id is unstable, then compute the simulated
%                      measurements matrix by uisng the one step ahead 
%                      prediction model 
%
%
% Meas_mtx_sim:  Simulated data matrix, generated by the identified model
%
% Err_mtx:  prediction error matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin~=4
    % Wrong number of input arguments
    error('Wrong number of input arguments');
end

[n_time_pts,n_nodes]=size(Meas_mtx);

num_exp=length(num_tps_4exp);
% num_tps_4exp(i)= number of time points for the i-th time-series

Meas_mtx_exp_cell{1,num_exp}=[];
Meas_mtx_sim_pr_exp_cell{1,num_exp}=[];
for num_ts=1:num_exp
    if num_ts>1
        ind_tps_in=sum(num_tps_4exp(1:num_ts-1))+1;
    else
        ind_tps_in=1;
    end
    
    ind_tps_end=sum(num_tps_4exp(1:num_ts));
      
    Meas_mtx_exp_cell{1,num_ts}=Meas_mtx(ind_tps_in:ind_tps_end,:);
    Meas_mtx_sim_pr_exp_cell{1,num_ts}=Meas_mtx_sim_pr(ind_tps_in:ind_tps_end,:);

 
end

IV_mtx_exp_cell{1,num_exp}=[];
Omega_mtx_exp_cell{1,num_exp}=[];
Theta_mtx_exp_cell{1,num_exp}=[];

for num_ts=1:num_exp
    
    Meas_mtx_sim_pr_exp= Meas_mtx_sim_pr_exp_cell{1,num_ts};
    Meas_mtx_exp= Meas_mtx_exp_cell{1,num_ts};
   
    IV_mtx_exp_cell{1,num_ts}= Meas_mtx_sim_pr_exp(1:end-1,:);
    Omega_mtx_exp_cell{1,num_ts}= Meas_mtx_exp(1:end-1,:);
    Theta_mtx_exp_cell{1,num_ts}= Meas_mtx_exp(2:end,:);

end


% Matrix of the measurements at time point k from 1 to n_time_pts-1
Omega_mtx=[];
% the IV matrix 
IV_mtx=[];
% Matrix of the measurements at time point k+1 from 2 to n_time_pts
Theta_mtx=[];
for num_ts=1:num_exp
    Omega_mtx_old=Omega_mtx;
    Omega_mtx_exp=Omega_mtx_exp_cell{1,num_ts};
    Omega_mtx=[Omega_mtx_old;Omega_mtx_exp];
    IV_mtx_old=IV_mtx;
    IV_mtx_exp=IV_mtx_exp_cell{1,num_ts};
    IV_mtx=[IV_mtx_old;IV_mtx_exp];
    Theta_mtx_old=Theta_mtx;
    Theta_mtx_exp=Theta_mtx_exp_cell{1,num_ts};
    Theta_mtx=[Theta_mtx_old;Theta_mtx_exp];
end

Omega_mtx=Omega_mtx';
Theta_mtx=Theta_mtx';
IV_mtx=IV_mtx';


U_mtx=zeros(num_exp,n_time_pts-num_exp);
for num_ts=1:num_exp
     if num_ts>1
        ind_tps_in=sum(num_tps_4exp(1:num_ts-1))-(num_ts-1)+1;
    else
        ind_tps_in=1;
    end
    
    ind_tps_end=sum(num_tps_4exp(1:num_ts))-num_ts;
    U_mtx(num_ts,ind_tps_in:ind_tps_end)=1;
end    
Omega_mtx=[Omega_mtx; U_mtx];
IV_mtx=[IV_mtx;U_mtx];
Rec_mtx=zeros(n_nodes,n_nodes+num_exp);

n_nodes_th_big=300;
% identification for each row beta=(X'*X+n*lambda*I)^(-1)*X'*y
% if use the IV method, X_tr=X'=IV_mtx; 
for idx=1:n_nodes
    
    X=Omega_mtx';
    y=Theta_mtx(idx,:)';
    % if Meas_mtx_sim_pr equals Meas_mtx then X_tr=X'                                     
    X_tr=IV_mtx;
    
    if not(isempty(null_constr_list))
        %find the elements that have to be nullified for each row
        idx_col2del=find(null_constr_list(:,1)==idx);

        if not(isempty(idx_col2del))
            %define a new linear regression, setting to zero the j-th
            %column of X
            col2del=null_constr_list(idx_col2del,2);
            X(:,col2del)=0;
            X_tr(col2del,:)=0;
        end
    end
    
    % define the Hessian
    H=X_tr*X;
    if n_nodes<n_nodes_th_big
        max_eig_H=max(abs(eig(H)));
        % choose a good ridge parameter
        %     lambda_min=fminbnd(@(x_min)lambda_GCV(x_min,X,H,X_tr,y),0,max_eig_H/size(X,1));
        
        lambda_min=fminbnd(@(x_min)lambda_GCV(x_min,X,H,X_tr,y),0,max_eig_H);
           
   
    else 
        if n_nodes==1565 % Ecoli network
            
            num_gen_pert=40; % number of genes perturbed for each perturbations
            
        else
            num_gen_pert=input('number of genes perturbed for each perturbation: '); 
        end
        if rem(idx-1,num_gen_pert)==0   % compute the  ridge parameter for each perturbation
            max_eig_H=max(abs(eig(H)));
            % choose a good ridge parameter
            %     lambda_min=fminbnd(@(x_min)lambda_GCV(x_min,X,H,X_tr,y),0,max_eig_H/size(X,1));
            
            lambda_min=fminbnd(@(x_min)lambda_GCV(x_min,X,H,X_tr,y),0,max_eig_H);
        end
        
    end
    beta=(H+size(X,1)*lambda_min*eye(size(H,1)))\X_tr*y;
    %     beta_ls=H\X'*y;
    %     Rec_mtx_ls(idx,:)=beta_ls';
    Rec_mtx(idx,:)=beta';
  
end


A_id=Rec_mtx(1:n_nodes,1:n_nodes);
B_id=Rec_mtx(1:n_nodes,n_nodes+1:n_nodes+num_exp);

% Compute the prediction error 
Xsim=Rec_mtx*Omega_mtx;

% Compute the error matrix
Err_mtx=Theta_mtx-Xsim;


% Compute the state evolution of the identified model and compare
% with experimental measurements

eig_Rec_mtx_LS=eig(A_id);
mod_eig_val=abs(eig_Rec_mtx_LS);
if max(mod_eig_val)>=1
    sol_found=false;
else
    sol_found=true;
end


% compute the "dynamic" error

Meas_mtx_sim=[];
Meas_mtx_pr=[];
norm_est_err_dy=zeros(1,num_exp);

for num_ts=1:num_exp
    %     Meas_mtx_r=zeros(n_time_pts,n_nodes);
    %     eval(['Meas_mtx_r=Meas_mtx_',num2str(num_ts),';']);
    n_time_pts_exp=num_tps_4exp(num_ts);
       
    if num_ts>1
        ind_tps_in=sum(num_tps_4exp(1:num_ts-1))-(num_ts-1)+1;
    else
        ind_tps_in=1;
    end
    
    ind_tps_end=sum(num_tps_4exp(1:num_ts))-num_ts;
    Meas_mtx_r= Meas_mtx_exp_cell{1,num_ts};
    [~,Meas_sim_ts]=dlsim(A_id,B_id(:,num_ts),zeros(1,n_nodes),0,ones(n_time_pts_exp,1),Meas_mtx_r(1,1:n_nodes)');
    
    Meas_mtx_pr_ts=[Meas_mtx_r(1,1:n_nodes);Xsim(:,ind_tps_in:ind_tps_end)'];
    Meas_mtx_pr_old=Meas_mtx_pr;
    Meas_mtx_pr=[Meas_mtx_pr_old; Meas_mtx_pr_ts];
    
    Meas_mtx_old=Meas_mtx_sim;
    Meas_mtx_sim=[Meas_mtx_old;Meas_sim_ts];
    Err_mtx_dy=abs(Meas_mtx_r-Meas_sim_ts);

    % Normalization of the error matrix
    % Find the maximum absolute value in the matrix
    max_Meas_mtx=max(max(abs(Meas_mtx_r)));
    % Divide each element by the abs value of the corresponding experimental
    % measurement
    for idx2=1:size(Meas_mtx_r,2)
        % The normalized error is computed only on the values that are at
        % least 10% of the maximum of the signal and only on those signals
        % whose maximum is at least 10% of the maximum over all the signals
        if max(abs(Meas_mtx_r(:,idx2)))>max_Meas_mtx*0.1
            col_thr=0.1*max(abs(Meas_mtx_r(:,idx2)));
        else
            % Do not consider this column in the computation of the estimation
            % error
            col_thr=Inf;
        end
        for idx1=1:size(Meas_mtx_r,1)
            if abs(Meas_mtx_r(idx1,idx2))>col_thr
                Err_mtx_dy(idx1,idx2)=Err_mtx_dy(idx1,idx2)/abs(Meas_mtx_r(idx1,idx2));
            else
                Err_mtx_dy(idx1,idx2)=0;
            end
        end
    end
    % Compute the cumulative error for each signal, pick the largest one and
    % divide by the number of time points
    norm_est_err_dy(num_ts)=norm(Err_mtx_dy,1)/n_time_pts_exp;
end

if not(sol_found)
    Meas_mtx_sim=Meas_mtx_pr;
end
