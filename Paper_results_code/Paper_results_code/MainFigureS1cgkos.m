%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script enables to compute the reconstruced connectivity matrices A,
% the perturbation vectors B and performance indexes obtained by Riva and 
% the other LS-based approaches for the 5 in-silico networks of 100 nodes 
% of Dream4 challenge, using the second part of the data-set, 
% when the perturbation is removed - 10 time-series, 
% each one of 11 time points; we replicate the data-set 20 times by adding 
% noise.   
% In the code, there is a flag (do_sim), which is true if you choose to do 
% the simulation, otherwise you can load the file 
% res_perf_B0_tps11_from_500to1000_rep20.mat 
% (within the folder Data_results_figs_3_4_S1/results/), where the 
% identified variables have been saved.
% The script generates the plots reported in the third column of Fig S1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
% 
% 
n_nodes=100; % test over the networks of twenty nodes 

num_exp=10; % number of time-series for exp.

% networks number
n_networks=5;

% number of tests
n_tests=21;




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

perf_aupr_LS=zeros(n_networks,n_tests);               
perf_auroc_LS=zeros(n_networks,n_tests);

perf_aupr_LS_IV=zeros(n_networks,n_tests);               
perf_auroc_LS_IV=zeros(n_networks,n_tests);

perf_aupr_RLS=zeros(n_networks,n_tests);               
perf_auroc_RLS=zeros(n_networks,n_tests);

perf_aupr_RLS_IV=zeros(n_networks,n_tests);               
perf_auroc_RLS_IV=zeros(n_networks,n_tests);

% % define a flag (do_sim), true if you choose to do the simulation,
% % otherwise you can load the file.mat, where the identified variables have 
% % been saved
do_sim=false;
% do_sim=true;

%                                     
if not(do_sim)
    
    
   load Data_results_figs_3_4_S1/results/res_perf_B0_tps11_from_500to1000_rep20.mat
  
    

else
    
    
    load  Data_results_figs_3_4_S1/Data/Dream4_100nodes_replica.mat
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%% COMPUTE THE PERFORMANCE OF LS-based techniques: LS,LS-IV, RLS, RIVA (RLS-IV) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for num_net=1:n_networks 

        % load the data for the current network

        eval(['load  Data_results_figs_3_4_S1/Data/DREAM4_GoldStandard_InSilico_Size100_',num2str(num_net),'.mat'])

        n_elem_not_null=nnz(original_network);
        spar_coef=(n_nodes*(n_nodes-1)-n_elem_not_null)/(n_nodes^2-n_nodes);
        %use the following formula if diag is not equal to zero(n_nodes^2-n_elem_not_null)/(n_nodes^2-n_nodes);
        spars_networks(num_net)=spar_coef; 

        for num_test=1:n_tests
            
            
            if num_test==1
            
          
                %v_tps=[1:11,22:32,43:53,64:74,85:95,106:116, 127:137, 148:158, 169:179, 190:200];
                 
                v_tps=[11:21,32:42,53:63,74:84,95:105,116:126, 137:147, 158:168, 179:189, 200:210];

            else
                
                
                %v_tps=[1:2:21,42:2:62,83:2:103,124:2:144,165:2:185,206:2:226, 247:2:267, 288:2:308, 329:2:349, 370:2:390];
              
                v_tps=[21:2:41,62:2:82,103:2:123,144:2:164,185:2:205,226:2:246, 267:2:287, 308:2:328, 349:2:369, 390:2:410];

        
               
            end
            Meas_mtx_all=rete{num_test,num_net};

            Meas_mtx=Meas_mtx_all(v_tps,2:end);
            B_mtx=zeros(n_nodes,num_exp);
            
            sign_constr_list=[];
        

            num_pts=11;
            num_tps_4exp=num_pts*ones(1,num_exp);
            
            % Compute the variables to be identified (the connectivity
            % matrix A and the exogenous perturbation matrix B)
            
            % Inference by using RLS and RIVA (RLS_IV)
            
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
            % matrix A and the exogenous perturbation matrix B)
                   
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
        
        save ('Data_results_figs_3_4_S1/results/res_perf_B0_tps11_from_500to1000_rep20_redo.mat')

        
        
        
    end
end



%%%%%% Plot the performance of RIVA and the other LS algorithms 

tag_riva='ob';
tag_genie='sk';
tag_bingo='dr';
tag_ls='^c';
tag_ls_iv='vm';
tag_rls='<g';



for num_id_net=1:n_networks
   
    
    
    
    aupr_Riva=perf_aupr_RLS_IV(num_id_net,:);
    auroc_Riva=perf_auroc_RLS_IV(num_id_net,:);
    
    aupr_LS=perf_aupr_LS(num_id_net,:);
    auroc_LS=perf_auroc_LS(num_id_net,:);
    
    aupr_LS_IV=perf_aupr_LS_IV(num_id_net,:);
    auroc_LS_IV=perf_auroc_LS_IV(num_id_net,:);
    
    aupr_RLS=perf_aupr_RLS(num_id_net,:);
    auroc_RLS=perf_auroc_RLS(num_id_net,:);
    
    f=figure;
    hold on 
    grid on
    plot(auroc_Riva,aupr_Riva,tag_riva,'linewidth',1.5);
    plot(auroc_LS,aupr_LS,tag_ls,'linewidth',1.5);
    plot(auroc_LS_IV,aupr_LS_IV,tag_ls_iv,'linewidth',1.5);

    plot(auroc_RLS,aupr_RLS,tag_rls,'linewidth',1.5);

    xlim([0.5 0.8])
    
    ylim([0 0.2])
    if num_id_net==1
        name_net='net. 1';
    elseif num_id_net==2
        name_net='net. 2';
    elseif num_id_net==3
        name_net='net. 3';
    elseif num_id_net==4
        name_net='net. 4';
    elseif num_id_net==5
        name_net='net. 5';
    end
    ylabel('AUPR')
    xlabel('AUROC')
    str_title=[name_net,', tps=',num2str(num_pts),' in [500 1000] min'];
    title(str_title)
    
%     title(name_net)
    set(f,'Position',[10 10 250 250])
    set(gca,'fontsize',12)
%     if num_id_net==5
%         
%         legend('RIVA','LS','LS-IV','RLS')
%     end
end

