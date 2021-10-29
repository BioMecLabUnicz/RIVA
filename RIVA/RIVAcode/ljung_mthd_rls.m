
function Meas_mtx_flt_all=ljung_mthd_rls(Meas_sim,Err_mtx,num_tps_4exp,ord_ar_no_model)

% filter the data (Meas_sim) with the autoregressive (AR) model, whose 
% order is defined by the scalar input (named ord_ar_no_model), for the 
% residuals obtained by the error matrix (Err_mtx). The input
% num_tps_4exp is a vector containing the number of time points for each 
% experiment (time-series), num_tps_4exp(i)= number of time points for the
% i-th time-series; the number of rows of Meas_mtx is equal to 
% sum(num_tps_4exp). The estimation of the coefficients of the AR filter 
% is performed through RLS.

[~,n_nodes]=size(Meas_sim);


% compute the residual from the error matrix
err_iv=[];
for idx=1:size(Err_mtx,2)
    err_iv_old=err_iv;
    err_iv=[err_iv_old Err_mtx(:,idx)'];
end

% estimate an AR  model for the residuals (err_iv) to extract
% the remaining information from err_iv:
% X_err = Z_err*par_ar_no_model + noise where X_err and Z_err are obtained 
% by the vector err_iv as follows

% size of each Err_mtx_num_ts for each time series
%num_tps_exp_i=num_tps_4exp(i)-1;
num_tps_min=min(num_tps_4exp);
if ord_ar_no_model>num_tps_min-1
    ord_ar_no_model=num_tps_min-1;
end

num_exp=length(num_tps_4exp);
Err_mtx_exp_cell{1,num_exp}=[];
err_iv_exp_cell{1,num_exp}=[];
for num_ts=1:num_exp
       
    if num_ts>1
        ind_tps_in=sum(num_tps_4exp(1:num_ts-1))-(num_ts-1)+1;
        ind_tps_er_in=(sum(num_tps_4exp(1:num_ts-1))-(num_ts-1))*n_nodes+1;
    else
        ind_tps_in=1;
         ind_tps_er_in=1;
    end
    
    ind_tps_end=sum(num_tps_4exp(1:num_ts))-num_ts;
    ind_tps_er_end=(sum(num_tps_4exp(1:num_ts))-num_ts)*n_nodes;
    Err_mtx_exp_cell{1,num_ts}=Err_mtx(:,ind_tps_in:ind_tps_end);
    err_iv_exp_cell{1,num_ts}=err_iv( ind_tps_er_in: ind_tps_er_end);
end

% vector that contains the columns from 2 to n_time_pts of each
% Err_mtx_num_ts
X_err=zeros(1,n_nodes*(sum(num_tps_4exp)-2*num_exp));
for num_ts=1:num_exp
    err_iv_exp=err_iv_exp_cell{1,num_ts};
    if num_ts>1
        ind_tps_x_in=(sum(num_tps_4exp(1:num_ts-1))-2*(num_ts-1))*n_nodes+1;
    else     
        ind_tps_x_in=1;
    end
    
    ind_tps_x_end=(sum(num_tps_4exp(1:num_ts))-2*num_ts)*n_nodes;
    X_err(ind_tps_x_in:ind_tps_x_end)=err_iv_exp(n_nodes+1:n_nodes*(num_tps_4exp(num_ts)-1));
    %     eval(['X_err((num_tps-1)*n_nodes*(num_ts-1)+1:(num_tps-1)*n_nodes*num_ts)= err_iv_',num2str(num_ts),...
    %         '(n_nodes+1:num_tps*n_nodes);'])
    
end
X_err=X_err';
Z_err=[];
for num_ts=1:num_exp
     
    num_tps=num_tps_4exp(num_ts)-1;
    Z_err_one=zeros((num_tps-1)*n_nodes,ord_ar_no_model);
    err_iv_exp=err_iv_exp_cell{1,num_ts};
   
    for idx=1:ord_ar_no_model
        
        Z_err_one(:,idx)=[zeros(1,(idx-1)*n_nodes) err_iv_exp(1:n_nodes*(num_tps-idx))];
    end
    Z_err_old=Z_err;
    Z_err=[Z_err_old;Z_err_one];
    clear Z_err_old
end


% define the Hessian
H=Z_err'*Z_err;

max_eig_H=max(eig(H));

    
if n_nodes<50
    
    lambda_min=fminbnd(@(x_min)lambda_GCV(x_min,Z_err,H,Z_err',X_err),0,max_eig_H);
    lambda_min_t=lambda_min*size(Z_err,1);

        
elseif n_nodes>=50 && n_nodes<300
    
    % choose a ridge parameter
    % for big networks use just few time series (num_exp_4_lam)
    %     num_exp_4_lam=2;
    %     ind_tps_x_end=(sum(num_tps_4exp(1:num_ts))-2*num_ts)*n_nodes
    %     Z_err_lam=Z_err(1:ind_tps_x_end,:);
    %     X_err_lam=X_err(1:ind_tps_x_end,:);
    %     H_lam=Z_err_lam'*Z_err_lam;
    %
    %     max_eig_H_lam=max(eig(H_lam));
    %     lambda_min=fminbnd(@(x_min)lambda_GCV(x_min,Z_err_lam,H_lam,Z_err_lam',X_err_lam),0,max_eig_H_lam/size(Z_err_lam,1));
    %     lambda_min_t=lambda_min*size(Z_err_lam,1);
    
    vct_lambda=zeros(1,num_exp);
    for num_ts=1:num_exp
        if num_ts>1
            ind_tps_x_in=(sum(num_tps_4exp(1:num_ts-1))-2*(num_ts-1))*n_nodes+1;
        else
            ind_tps_x_in=1;
        end
        
        ind_tps_x_end=(sum(num_tps_4exp(1:num_ts))-2*num_ts)*n_nodes;
    
        Z_err_lam=Z_err(ind_tps_x_in:ind_tps_x_end,:);
        X_err_lam=X_err(ind_tps_x_in:ind_tps_x_end,:);
        H_lam=Z_err_lam'*Z_err_lam;

        max_eig_H_lam=max(eig(H_lam));
        lambda_min=fminbnd(@(x_min)lambda_GCV(x_min,Z_err_lam,H_lam,Z_err_lam',X_err_lam),0,max_eig_H_lam);
        lambda_min_exp=lambda_min*size(Z_err_lam,1);
        vct_lambda(num_ts)=lambda_min_exp;
    end
    lambda_min_t=median(vct_lambda);
else
    
    max_eig_H=max(abs(eig(H)));
    min_eig_H=min(abs(eig(H)));
    
    lambda_min_t=compute_lambada(min_eig_H,max_eig_H);  
    
end



% estimate the filter parameters
par_ar_no_model=(H+lambda_min_t*eye(size(H,1)))\Z_err'*X_err;

par_ar_no_model=-par_ar_no_model;

Meas_mtx_flt_all=[];
for num_ts=1:num_exp
    
    n_time_pts_4exp=num_tps_4exp(num_ts);
    % Filter the simulated Meas_sim with the AR filter
    Meas_sim_flt=zeros(n_time_pts_4exp,n_nodes);
     
    if num_ts>1
        ind_tps_in=sum(num_tps_4exp(1:num_ts-1))+1;
    else
        ind_tps_in=1;
    end
    ind_tps_end=sum(num_tps_4exp(1:num_ts));
    Meas_sim_ts=Meas_sim(ind_tps_in:ind_tps_end,:);
    for idx1=1:n_nodes
        for idx2=1:n_time_pts_4exp
            if idx2==1
                Meas_sim_flt(idx2,idx1)=Meas_sim_ts(idx2,idx1);
            elseif idx2>=2 && idx2<=ord_ar_no_model
                Meas_sim_flt(idx2,idx1)=Meas_sim_ts(idx2,idx1);
                for idx3=1:idx2-1
                    Meas_sim_flt_prev=Meas_sim_flt(idx2,idx1);
                    Meas_sim_flt(idx2,idx1)=Meas_sim_flt_prev+par_ar_no_model(idx3)*Meas_sim_ts(idx2-idx3,idx1);
                end
            elseif idx2>ord_ar_no_model
                Meas_sim_flt(idx2,idx1)=Meas_sim_ts(idx2,idx1);
                for idx3=1:ord_ar_no_model
                    Meas_sim_flt_prev=Meas_sim_flt(idx2,idx1);
                    Meas_sim_flt(idx2,idx1)=Meas_sim_flt_prev+par_ar_no_model(idx3)*Meas_sim_ts(idx2-idx3,idx1);
                end
            end
        end
    end
    Meas_mtx_flt_all_old=Meas_mtx_flt_all;
    Meas_mtx_flt_all=[Meas_mtx_flt_all_old;Meas_sim_flt];
end

