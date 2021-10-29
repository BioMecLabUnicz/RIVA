%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script enables to compute the reconstruced connectivity matrices A,
% the perturbation matrices B and performance indexes obtained by Riva for 
% the 5 in-silico networks of 100 nodes provided in the DREAM3 Challenge,
% which are used for the generation of our dataset through GeneNetWeaver. 
% The script generates the results shown in Figure 1 of the ms. 
% In the code, there is a flag (do_sim), which is true if you choose to do 
% the simulation, otherwise you can load the file 
% res_perf_Bempty_InSilicoSize100_tps21_rep5__multifactorial_timeseries.mat
% (within the folder Data_results_fig1/results), where the identified 
% variables have been saved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
% 
% 
n_nodes=100; % test over the networks of 100 nodes 

num_exp=10; % number of time-series for experiment

type_data='_multifactorial_timeseries';

% networks number
n_networks=5;

% number of tests
n_tests=5;

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
%  do_sim=true;

%                                     
if not(do_sim)
    
    
   eval(['load Data_results_fig1/results/res_perf_Bempty_InSilicoSize100_tps21_rep5_',type_data])
  

else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%% COMPUTE THE PERFORMANCE OF LS-based techniques: LS, LS-IV, RLS, RIVA (RLS-IV) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for num_net=1:n_networks 


        if num_net==1
            name_net='Ecoli1';
        elseif num_net==2
            name_net='Ecoli2';
        elseif num_net==3
            name_net='Yeast1';
        elseif num_net==4
            %n_test_in=17;
            name_net='Yeast2';
        elseif num_net==5
            name_net='Yeast3';
            %n_test_in=1;
        end

        % load the data for the current network
    % 
        eval(['load Data_results_fig1/',name_net,'/InSilicoSize100-',name_net,type_data,'.mat']);

        eval(['load Data_results_fig1/',name_net,'/InSilicoSize100-',name_net,'_goldstandard.mat']);


        n_elem_not_null=nnz(original_network);
        spar_coef=(n_nodes*(n_nodes-1)-n_elem_not_null)/(n_nodes^2-n_nodes);
        %use the following formula if diag is not equal to zero(n_nodes^2-n_elem_not_null)/(n_nodes^2-n_nodes);
        spars_networks(num_net)=spar_coef; 


        for num_test=1:n_tests

            

            num_pts=21;
            eval(['Meas_mtx=Meas_mtx',type_data,'(num_pts*num_exp*(num_test-1)+1:num_pts*num_exp*num_test,:);']);
            B_mtx=[];
            sign_constr_list=[];
            num_tps_4exp=num_pts*ones(1,num_exp);
            
            % Compute the variables to be identified (the connectivity
            % matrix A and the exogenous perturbation matrix B)
            
            % Inference by using RLS and RIVA (RLS with IV)  
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
        
        
        eval(['save Data_results_fig1/results/res_perf_Bempty_InSilicoSize100_tps21_rep5_',type_data,'_redo'])
        
    end
end


%%%%%% Plot the performance of RIVA and the other algorithms


load Data_results_fig1/results/res_Bingo
load Data_results_fig1/results/res_dynGenie3

tag_riva='ob';
tag_genie='sk';
tag_bingo='dr';
tag_ls='^c';
tag_ls_iv='vm';
tag_rls='<g';


aupr_mtx_Riva=perf_aupr_RLS_IV;
auroc_mtx_Riva=perf_auroc_RLS_IV;
for num_id_net=1:n_networks
   
    
    if num_id_net==1
            name_net='Ecoli1';
        elseif num_id_net==2
            name_net='Ecoli2';
        elseif num_id_net==3
            name_net='Yeast1';
        elseif num_id_net==4
            name_net='Yeast2';
        elseif num_id_net==5
            name_net='Yeast3';
    end
    
    aupr_Riva=aupr_mtx_Riva(num_id_net,:);
    auroc_Riva=auroc_mtx_Riva(num_id_net,:);
   
    f=figure;
    hold on 
    grid on
    plot(auroc_Riva,aupr_Riva,tag_riva,'linewidth',1.5);
    plot(Results_dynGenie3_100.auroc(num_id_net,:),Results_dynGenie3_100.aupr(num_id_net,:),tag_genie,'linewidth',1.5);
    plot(Results_Bingo_100.auroc(num_id_net,:),Results_Bingo_100.aupr(num_id_net,:),tag_bingo,'linewidth',1.5);
    xlim([0.5 0.8])
     
    ylim([0 0.105])
    ylabel('AUPR')
    xlabel('AUROC')
    title(name_net)
    set(f,'Position',[10 10 250 250])
    set(gca,'fontsize',12)
    if num_id_net==5
        
        legend('RIVA','dynGenie3','Bingo')
    end
end


