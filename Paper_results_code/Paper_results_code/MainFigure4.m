
% This scrip allows to generate the results shown in Figure 4 for the five 
% 100-node networks of DREAM4 challenge.
% Compute the perfomance by evaluating the final list of edges on the
% basis of the edge average ranking over the replicated inference tests.
% We used only the first half of the time-series data provided by DREAM4, 
% without considering the half interval without perturbation.
% We replicated the dataset by adding noise, moreover we generated replicas 
% with an increasing number of points by interpolating the dataset.

% load the stored identified networks obtained by running
% MainFigureS1aeimq.m - i.e.identification using the first half
% of data - 10 time series of 11 time points - we replicate the data-set 20
% times by adding noise

load Data_results_figs_3_4_S1/results/res_perf_Bempty_tps11_from_0to500_rep20.mat

A_id_cell_RLS_IV_tot=A_id_cell_RLS_IV;


% load the stored identified networks obtained by running
% MainFigureS1bfjnr.m - i.e.identification using the first half
% of data - 10 time series of 11 time points - we replicate the data-set 20
% times by adding noise and sampling the data-set for generating time series 
% of 21 time points

load Data_results_figs_3_4_S1/results/res_perf_Bempty_tps21_from_0to500_rep20.mat

% connect the two set of reconstructed networks
A_id_cell_RLS_IV_tot=[A_id_cell_RLS_IV_tot A_id_cell_RLS_IV(:,2:end)];

aupr_Riva=zeros(1,n_networks);
auroc_Riva=zeros(1,n_networks);

n_networks=size(A_id_cell_RLS_IV_tot,1);
all_tests=size(A_id_cell_RLS_IV_tot,2);

% number of links to be added for computing performance at each evaluation (idx_eval): from 1 to 9990=sum(inc_coef2add); 
inc_coef2add=[1 2 2 3 3 4:35 40 50 75 100 200 500 800 1000 1500 2000 3000];
num_mode=1; %criterion of normalization for ranking, it can be also 2 or 3, see norm_rec_mtx.m


for num_id_net=1:n_networks
 
    Constr_list_cell_net=[];
    for num_test=1:all_tests
        
        A_id_cell=A_id_cell_RLS_IV_tot{num_id_net,num_test};
        %normalize the identified matrix
        Normal_A_id = norm_rec_mtx(A_id_cell,num_mode);
        for idx2=1:n_nodes
            Normal_A_id(idx2,idx2)=0;
        end
        % sort the matrix coef.
        [idx_row, idx_col, val_coef]=find(Normal_A_id);
        [val_coef_sort, idx_sort]=sort(val_coef,'descend');
        %create the ranking list
        Constr_list=[idx_row(idx_sort) idx_col(idx_sort) val_coef_sort];
        
        % store the ranking list for a test
        Constr_list_cell_net{num_test}=Constr_list;
        
        
    end
    

    
    %load the original network
    eval(['load  Data_results_figs_3_4_S1/Data/DREAM4_GoldStandard_InSilico_Size100_',num2str(num_id_net),'.mat'])

    num_eval=length(inc_coef2add);
    PPV=zeros(1,num_eval);
    Se=zeros(1,num_eval);
    Sp=zeros(1,num_eval);
    FP=zeros(1,num_eval);
    TN=(n_nodes-1)*n_nodes-nnz(original_network);

    for idx_eval=1:num_eval
        
        A_n_rank=zeros(n_nodes,n_nodes); % matrix used to create an average ranking of each edge over the inference tests
        num_edges=sum(inc_coef2add(1:idx_eval));% number of edges to be inserted                  
        if not(isempty(Constr_list_cell_net))            
            
            for idx_id=1:size(Constr_list_cell_net,2) % for all the tests
                rank_list=Constr_list_cell_net{idx_id}; 
               
                % for the entry of A_n_rank add 1 if the corresponding
                % edge is in the ranking list (of n edges, n=num_edges) for
                % the given test
                for idx_n=1:num_edges
                    
                    A_n_rank(rank_list(idx_n,1),rank_list(idx_n,2))=A_n_rank(rank_list(idx_n,1),rank_list(idx_n,2))+1;
                    
                    
                end
                
                
            end
        end
        
        A_n_rank_cell{num_id_net,idx_eval}=A_n_rank;
        A_id_ones=val_rec_list(A_n_rank,num_edges); % A matrix with entries equal to 1 corresponding to the n edges (n=num_edges) with high probability 
        [Perf_idx, Res]= ppv_sens_spec(original_network,A_id_ones);
        PPV(idx_eval)=Perf_idx(2,1);
        Se(idx_eval)=Perf_idx(2,2);
        Sp(idx_eval)=Perf_idx(2,3);
        FP(idx_eval)=size(Res.Dir.False_pos,1);
        
    end
            
    aupr_Riva(num_id_net)=eval_aupr(PPV,Se);
    auroc_Riva(num_id_net)=eval_auroc(Se,FP,TN);
                       
end


%%%%%% Plot the performance of RIVA and the other algorithms 


load Data_results_figs_3_4_S1/results/Bingo_100
load Data_results_figs_3_4_S1/results/dynGenie3_100

tag_riva='ob';
tag_genie='sk';
tag_bingo='dr';

for num_id_net=1:5
   
    
   
    
    if num_id_net==1
            name_net='net. 1';
        elseif num_id_net==2
            name_net='net. 2';
        elseif num_id_net==3
            name_net='net. 3';
        elseif num_id_net==4
            %n_test_in=17;
            name_net='net. 4';
        elseif num_id_net==5
            name_net='net. 5';
            %n_test_in=1;
    end
%     
    
    f=figure;
    hold on 
    grid on
    plot(auroc_Riva(num_id_net),aupr_Riva(num_id_net),tag_riva,'linewidth',1.5);
    plot(Result_dynGenie3100.auroc(num_id_net),Result_dynGenie3100.aupr(num_id_net),tag_genie,'linewidth',1.5);
    plot(Result_Bingo100.auroc(num_id_net),Result_Bingo100.aupr(num_id_net),tag_bingo,'linewidth',1.5);
    xlim([0.5 0.8])
     
    ylim([0 0.15])
    ylabel('AUPR')
    xlabel('AUROC')
    title(name_net)
    set(f,'Position',[10 10 250 250])
    set(gca,'fontsize',12)
    
     if num_id_net==5
        
        legend('RIVA','dynGenie3','Bingo')
    end
    
end



