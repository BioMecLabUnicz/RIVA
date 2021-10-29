%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script enables to compute the reconstruced cell cycle regulatory 
% subnetwork in Saccharomyces cerevisiae  from experimental microarray data
% and evaluate the performance obtained by RIVA and compared with Bingo and
% dynGenie3 (Fig. 5 and 6, cases a) b) and c)).
% In the code, there is a flag (do_sim), which is true if you choose to do 
% the simulation, otherwise you can load the file where the identified 
% variables have been saved:
% emtab_643_15_res.mat, emtab_643_10_res.mat or 
% emtab_643_8_res.mat 
% (within the folder Bio_Example_data_results/results/), depending on
% the selected sample time ts (within the script, set the parameter 
% val_ts= 15, 10 or 8).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
 clc
 clear all
 close all
 
% folder name containing data and results 
folder_exp='Bio_Example_data_results';

% % load a cell array containig the standard (first column) and systematic 
% % (second column) names of genes that compose the 
% % subnetwork
 load Bio_Example_data_results/name_genes_subnetwork
% 
% load the data 
% the sample time is 15, 10 or 8, set val_ts=15, 10 or 8
val_ts=15; 

file_name_1=['emtab_643_',num2str(val_ts)];
   
%file_name_2=['emtab_1908_L_',num2str(val_ts)];

 
eval (['load Bio_Example_data_results/',file_name_1,'.mat'])

% eval (['load Bio_Example_data_results/',file_name_2,'.mat'])


% % load the original network to identify, obtained by the BioGRID database, 
% % taking into account the information of Simon et al  (2001) and Bahler 
% % (2005). 
 load Bio_Example_data_results/true_subnetwork

% define the number of nodes of the subnetwork
n_nodes=size(original_network,1);


% define the sparsity coefficient for the subnetwork
spar_coef=(n_nodes^2-(nnz(original_network(1:27,1:27))-...
    nnz(diag(original_network(1:27,1:27)))))/(n_nodes^2-n_nodes);



% number id used for store the results 
n_networks=1;

% number of tests
n_tests=50; % for computing performance we use the first 20 tests

num_mode=1; % type of normalization for ranking, see norm_rec_mtx.m


% allocate a memory space to save the output variables of the iterative
% procedure

A_id_cell_RLS{n_networks,n_tests}=[]; % cell array containing the identified 
                                     % connectivity matrices, obtained by
                                     % RLS, for each network and 
                                     % test  

A_id_cell_RLS_IV{n_networks,n_tests}=[]; % cell array containing the identified 
                                     % connectivity matrices, obtained by
                                     % RIVA (RLS with IV), for each network and 
                                     % test  

                                     
B_id_cell_RLS{n_networks,n_tests}=[]; % cell array containing the identified 
                                     % input matrices, obtained by
                                     % RLS, for each network and 
                                     % test  

B_id_cell_RLS_IV{n_networks,n_tests}=[]; % cell array containing the identified 
                                     % input matrices, obtained by
                                     % RIVA, for each network and test
                                       
                                     
est_err_cell_RLS{n_networks,n_tests}=[]; % cell array containing the errors 
                                         % of the identified models, obtained by
                                     % RLS, for each network and test  

est_err_cell_RLS_IV{n_networks,n_tests}=[]; % cell array containing the errors 
                                         % of the identified models, obtained by
                                     % RIVA, for each network and test
                                                                          
sol_found_RLS{n_networks,n_tests}=[]; % cell array containing a flag on the 
                                     % stability of the identified systems, determining 
                                     % the mode for obtaining the simulated data matrices   
                                     % by RLS, for each network and test

sol_found_RLS_IV{n_networks,n_tests}=[]; % cell array containing a flag on the 
                                     % stability of the identified systems, determining 
                                     % the mode for obtaining the simulated data matrices   
                                     % by RIVA, for each network and test
                                     
                                                                          
A_id_cell_LS{n_networks,n_tests}=[]; % cell array containing the identified 
                                     % connectivity matrices, obtained by
                                     % LS, for each network and test  

A_id_cell_LS_IV{n_networks,n_tests}=[]; % cell array containing the identified 
                                     % connectivity matrices, obtained by
                                     % LS-IV, for each network and test
                                  
B_id_cell_LS{n_networks,n_tests}=[]; % cell array containing the identified 
                                     % input matrices, obtained by
                                     % LS, for each network and test
                                  
B_id_cell_LS_IV{n_networks,n_tests}=[]; % cell array containing the identified 
                                     % input matrices, obtained by
                                     % LS-IV, for each network and test
                               
est_err_cell_LS{n_networks,n_tests}=[]; % cell array containing the errors 
                                         % of the identified models, obtained by
                                     % LS, for each network and test
                                
est_err_cell_LS_IV{n_networks,n_tests}=[]; % cell array containing the errors 
                                         % of the identified models, obtained by
                                     % LS-IV, for each network and test
                                     
sol_found_LS{n_networks,n_tests}=[]; % cell array containing a flag on the 
                                     % stability of the identified systems, determining 
                                     % the mode for obtaining the simulated data matrices   
                                     % by LS, for each network and test

sol_found_LS_IV{n_networks,n_tests}=[]; % cell array containing a flag on the 
                                     % stability of the identified systems, determining  
                                     % the mode for obtaining the simulated data matrices   
                                     % by LS-IV, for each network and test 

% define a vector to store the sparsity coefficient for each network
spars_networks=zeros(1,n_networks);

% define array to store the results 
PPV_LS{n_networks,n_tests}=[];
Se_LS{n_networks,n_tests}=[];
Sp_LS{n_networks,n_tests}=[];
FP_LS{n_networks,n_tests}=[];        

PPV_RLS{n_networks,n_tests}=[];
Se_RLS{n_networks,n_tests}=[];
Sp_RLS{n_networks,n_tests}=[];
FP_RLS{n_networks,n_tests}=[];        


PPV_RLS_IV{n_networks,n_tests}=[];
Se_RLS_IV{n_networks,n_tests}=[];
Sp_RLS_IV{n_networks,n_tests}=[];
FP_RLS_IV{n_networks,n_tests}=[];        

PPV_LS1{n_networks,n_tests}=[];
Se_LS1{n_networks,n_tests}=[];
Sp_LS1{n_networks,n_tests}=[];
FP_LS1{n_networks,n_tests}=[];        


PPV_LS1_IV{n_networks,n_tests}=[];
Se_LS1_IV{n_networks,n_tests}=[];
Sp_LS1_IV{n_networks,n_tests}=[];
FP_LS1_IV{n_networks,n_tests}=[];        

perf_aupr_LS=NaN(n_networks,n_tests);               
perf_auroc_LS=NaN(n_networks,n_tests);

perf_aupr_LS_IV=NaN(n_networks,n_tests);               
perf_auroc_LS_IV=NaN(n_networks,n_tests);

perf_aupr_RLS=NaN(n_networks,n_tests);               
perf_auroc_RLS=NaN(n_networks,n_tests);

perf_aupr_RLS_IV=NaN(n_networks,n_tests);               
perf_auroc_RLS_IV=NaN(n_networks,n_tests);

% % define a flag (do_sim), true if you choose to do the simulation,
% % otherwise you can load the file.mat, where the identified variables have 
% % been saved
%do_sim=true;
do_sim=false;
%                                     
if not(do_sim)
    
    filename= ['Bio_Example_data_results/results/',file_name_1,'_res.mat'];
    load  (filename)
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% COMPUTE THE PERFORMANCE OF LS-based techniques: LS,LS-IV, RLS, RIVA (RLS-IV) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for num_net=1
        
       
        B_mtx=[];
        sign_constr_list=[];
        

        str_array_2=['emtab_643_',num2str(val_ts)];

        
        %         str_array=['emtab_1908_L_',num2str(val_ts)];
        %         time_1=5;
        %         time_2=105;
        %         t_sample=val_ts;
        %         v_time=5:t_sample:205;
        % find the samples for the selected time interval, from 5 to 105 minutes.
        %         [idx_time1,idx_time2]=find_idx_time_4ts(time_1,time_2,v_time);

        for num_test=1:n_tests
                       
            
            eval (['Meas_mtx_2=',str_array_2,'{1,num_test};']);
            
            
            Meas_mtx_2=Meas_mtx_2';      
            
            
            Meas_mtx=Meas_mtx_2;
             
        
            % number of time points for experiment
            num_tps_4exp(1)=size(Meas_mtx,1);
            
                
             % Compute the variables to be identified (the connectivity
            % matrix A and the exogenous perturbation vector B)
            
            % Inference by using RLS and RLS_IV
            
            [A_id_cell,B_id_cell,est_err_cell,sol_found] = ...
                net_rec_riva(Meas_mtx,B_mtx,sign_constr_list,num_tps_4exp);
            
            
            % store the identified variables by RLS          
            A_id_cell_RLS{num_net,num_test}=A_id_cell{1,1};
            B_id_cell_RLS{num_net,num_test}=B_id_cell{1,1};
            est_err_cell_RLS{num_net,num_test}=est_err_cell{1,1};
            sol_found_RLS{num_net,num_test}=sol_found{1,1};
            
            % store the identified variables by RIVA
            A_id_cell_RLS_IV{num_net,num_test}=A_id_cell{1,2};
            B_id_cell_RLS_IV{num_net,num_test}=B_id_cell{1,2};
            est_err_cell_RLS_IV{num_net,num_test}=est_err_cell{1,2};
            sol_found_RLS_IV{num_net,num_test}=sol_found{1,2};
          

            
            % Compute the variables to be identified (the connectivity
            % matrix A and the exogenous perturbation vector B)
                   
            % Inference by using LS and LS-IV
            
            [A_id_cell,B_id_cell,est_err_cell,sol_found] = ...
                net_rec_ls_iv(Meas_mtx,B_mtx,sign_constr_list,num_tps_4exp);

            
           
            % store the identified variables by LS
            A_id_cell_LS{num_net,num_test}=A_id_cell{1,1};
            B_id_cell_LS{num_net,num_test}=B_id_cell{1,1};
            est_err_cell_LS{num_net,num_test}=est_err_cell{1,1};
            sol_found_LS{num_net,num_test}=sol_found{1,1};
            
            % store the identified variables by LS-IV
            A_id_cell_LS_IV{num_net,num_test}=A_id_cell{1,2};
            B_id_cell_LS_IV{num_net,num_test}=B_id_cell{1,2};
            est_err_cell_LS_IV{num_net,num_test}=est_err_cell{1,2};
            sol_found_LS_IV{num_net,num_test}=sol_found{1,2};
            




            % evaluation of the performance by calculating the AUPR,
            % AUROC, PPV, Sensitivity, Specificity and False Positives
            
            
            % by LS
            [AUPR_m, AUROC_m, PPV, Se, Sp, FP]=val_aupr_auroc(original_network,A_id_cell_LS{num_net,num_test},num_mode);
            
            % store LS performance
            perf_aupr_LS(num_net,num_test)=AUPR_m;
            perf_auroc_LS(num_net,num_test)=AUROC_m;
            
            PPV_LS{num_net,num_test}=PPV;
            Se_LS{num_net,num_test}=Se;
            Sp_LS{num_net,num_test}=Sp;
            FP_LS{num_net,num_test}=FP;
            
            clear AUPR_m AUROC_m PPV Se Sp FP
            
            % by LS-IV
            [AUPR_m, AUROC_m, PPV, Se, Sp, FP]=val_aupr_auroc(original_network,A_id_cell_LS_IV{num_net,num_test},num_mode);
            
            % store LS-IV performance
            perf_aupr_LS_IV(num_net,num_test)=AUPR_m;
            perf_auroc_LS_IV(num_net,num_test)=AUROC_m;
            
            PPV_LS1_IV{num_net,num_test}=PPV;
            Se_LS1_IV{num_net,num_test}=Se;
            Sp_LS1_IV{num_net,num_test}=Sp;
            FP_LS1_IV{num_net,num_test}=FP;
            
            clear AUPR_m AUROC_m PPV Se Sp FP
            
            % by RLS
            [AUPR_m, AUROC_m, PPV, Se, Sp, FP]=val_aupr_auroc(original_network,A_id_cell_RLS{num_net,num_test},num_mode);
            
            % store RLS performance
            perf_aupr_RLS(num_net,num_test)=AUPR_m;
            perf_auroc_RLS(num_net,num_test)=AUROC_m;
            
            PPV_RLS{num_net,num_test}=PPV;
            Se_RLS{num_net,num_test}=Se;
            Sp_RLS{num_net,num_test}=Sp;
            FP_RLS{num_net,num_test}=FP;
            
            clear AUPR_m AUROC_m PPV Se Sp FP
            
            % by RIVA
            [AUPR_m, AUROC_m, PPV, Se, Sp, FP]=val_aupr_auroc(original_network,A_id_cell_RLS_IV{num_net,num_test},num_mode);
            
            % store RIVA performance
            perf_aupr_RLS_IV(num_net,num_test)=AUPR_m;
            perf_auroc_RLS_IV(num_net,num_test)=AUROC_m;
            
            PPV_RLS_IV{num_net,num_test}=PPV;
            Se_RLS_IV{num_net,num_test}=Se;
            Sp_RLS_IV{num_net,num_test}=Sp;
            FP_RLS_IV{num_net,num_test}=FP;
            
            
            clear AUPR_m AUROC_m PPV Se Sp FP
            
            
          
        end
        
      
         eval (['save Bio_Example_data_results/results/',file_name_1,'_redo_res.mat'])
        
        
    end
    
end

%%%%%% Plot the performance of RIVA and the other algorithms 


%%%% tag for the experimental case (=1 for EMTAB 643, =2 for EMTAB 1908, =3
%%%% for both)
tag_case=1;

% Riva performance
AUROC_val=(perf_auroc_RLS_IV(1,:));
                        
AUPR_val=(perf_aupr_RLS_IV(1,:));

if val_ts==8
    
    load('Bio_Example_data_results/results/Bingo_8.mat')
    load ('Bio_Example_data_results/results/dynGenie3_8.mat')
    str_ts='8 min.';
    aupr_Bingo=Result_Bingo_8.aupr(tag_case,:);
    auroc_Bingo=Result_Bingo_8.auroc(tag_case,:);
    
    aupr_dynGenie3=Result_dynGenie3_8.aupr(tag_case,:);
    auroc_dynGenie3=Result_dynGenie3_8.auroc(tag_case,:);
    
    
    
elseif val_ts==10
    load('Bio_Example_data_results/results/Bingo_10.mat')
    load ('Bio_Example_data_results/results/dynGenie3_10.mat')
    str_ts='10 min.';
    aupr_Bingo=Result_Bingo_10.aupr(tag_case,:);
    auroc_Bingo=Result_Bingo_10.auroc(tag_case,:);
    
    aupr_dynGenie3=Result_dynGenie3_10.aupr(tag_case,:);
    auroc_dynGenie3=Result_dynGenie3_10.auroc(tag_case,:);
    
    
    
elseif val_ts==15
    load('Bio_Example_data_results/results/Bingo_15.mat')
    load ('Bio_Example_data_results/results/dynGenie3_15.mat')
    str_ts='15 min.';
    aupr_Bingo=Result_Bingo_15.aupr(tag_case,:);
    auroc_Bingo=Result_Bingo_15.auroc(tag_case,:);
    
    aupr_dynGenie3=Result_dynGenie3_15.aupr(tag_case,:);
    auroc_dynGenie3=Result_dynGenie3_15.auroc(tag_case,:);
    
end

if tag_case==1
    
    str_data='E-MTAB-643';
elseif tag_case==2
    str_data='E-MTAB-1908';
elseif tag_case==3
    str_data='E-MTAB-643:1908';
end
tag_col_mark='ob';
tag_genie='sk';
tag_bingo='dr';
f=figure;
hold on
grid on
str_title=[str_data, '; t_s=',str_ts];

aupr_Bingo=aupr_Bingo(1:20);
auroc_Bingo=auroc_Bingo(1:20);
aupr_dynGenie3=aupr_dynGenie3(1:20);
auroc_dynGenie3=auroc_dynGenie3(1:20);

AUPR_val=AUPR_val(1:20);
AUROC_val=AUROC_val(1:20);

plot((auroc_Bingo),(aupr_Bingo),tag_bingo);
plot((auroc_dynGenie3),(aupr_dynGenie3),tag_genie);
plot(AUROC_val,AUPR_val,tag_col_mark);

ylabel('AUPR')
xlabel('AUROC')
title(str_title)
    
set(f,'Position',[10 10 250 250])
  
xlim([0.4 .62])
ylim([0.13 0.25])


f=figure;

hold on
grid on

errorbar(mean(auroc_Bingo),mean(aupr_Bingo),std(aupr_Bingo),std(aupr_Bingo),std(auroc_Bingo),std(auroc_Bingo),tag_bingo,'linewidth',2);
errorbar(mean(auroc_dynGenie3),mean(aupr_dynGenie3),std(aupr_dynGenie3),std(aupr_dynGenie3),std(auroc_dynGenie3),std(auroc_dynGenie3),tag_genie,'linewidth',2);

errorbar(mean(AUROC_val),mean(AUPR_val),std(AUPR_val),std(AUPR_val),std(AUROC_val),std(AUROC_val),tag_col_mark,'linewidth',2);

ylabel('AUPR')
xlabel('AUROC')
title(str_title)
    
set(f,'Position',[10 10 250 250])
  
xlim([0.4 .62])
ylim([0.13 0.25])


