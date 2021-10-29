% This scrip allows to generate the results shown in Figure 3 for the five 
% 100-node networks of DREAM4 challenge.
% Compute the RIVA perfomance by evaluating the final list of edges on the
% basis of the edge average ranking over the replicated inference tests.
% We replicated the Dream4 challenge data (all the dataset) 
% where perturbations are applied on some genes/nodes of the network
% for the first half of the recording time and then are removed for the
% second half. We also generated replicas with an increasing number of 
% points by interpolating the dataset.


num_cases=2;  % we assume two different identifications by divinding the data-set by two: 
              % num_case=1, for data with perturbation from 0 to 500,
              % first half of data;
              % num_case=2 for data without perturbation from 500 to 1000,
              % second half of data.


% load the stored identified networks obtained by running
% MainFigureS1aeimq.m - i.e.identification using the first half
% of data - 10 time series of 11 time points - we replicate the data-set 20
% times by adding noise

load Data_results_figs_3_4_S1/results/res_perf_Bempty_tps11_from_0to500_rep20.mat

A_id_cell_RLS_IV_1=A_id_cell_RLS_IV;


% load the stored identified networks obtained by running
% MainFigureS1bfjnr.m - i.e.identification using the first half
% of data - 10 time series of 11 time points - we replicate the data-set 20
% times by adding noise and sampling the data-set for generating time series 
% of 21 time points

load Data_results_figs_3_4_S1/results/res_perf_Bempty_tps21_from_0to500_rep20.mat

% connect the two set of reconstructed networks
A_id_cell_RLS_IV_1=[A_id_cell_RLS_IV_1 A_id_cell_RLS_IV(:,2:end)];


% load the stored identified networks obtained by running
% MainFigureS1cgkos.m - i.e.identification using the second half
% of data - 10 time series of 11 time points - we replicate the data-set 20
% times by adding noise 

load Data_results_figs_3_4_S1/results/res_perf_B0_tps11_from_500to1000_rep20.mat
  
A_id_cell_RLS_IV_2= A_id_cell_RLS_IV;



% load the stored identified networks obtained by running
% FigureS1dhlpt.m - i.e.identification using the second half
% of data - 10 time series of 11 time points - we replicate the data-set 20
% times by adding noise and sampling the data-set for generating time series 
% of 21 time points

load Data_results_figs_3_4_S1/results/res_perf_B0_tps21_from_500to1000_rep20.mat
  
% connect the two set of reconstructed networks
A_id_cell_RLS_IV_2=[A_id_cell_RLS_IV_2 A_id_cell_RLS_IV(:,2:end)];


% 
aupr_Riva=zeros(1,n_networks);
auroc_Riva=zeros(1,n_networks);

n_networks=size(A_id_cell_RLS_IV_1,1);
all_tests=size(A_id_cell_RLS_IV_1,2);

% number of links to be added for computing performance at each evaluation (idx_eval): from 1 to 9990=sum(inc_coef2add); 
inc_coef2add=[1 2 2 3 3 4:35 40 50 75 100 200 500 800 1000 1500 2000 3000];
num_mode=1; %criterion of normalization for ranking, see norm_rec_mtx.m

                   
for num_id_net=1:n_networks
    
    Constr_list_cell_net=[];
    for num_case=1:num_cases
        
        
        if num_case==1
            
           
            A_id_cell_RLS_IV_case=A_id_cell_RLS_IV_1;
            
        elseif num_case==2
            
            A_id_cell_RLS_IV_case=A_id_cell_RLS_IV_2;
            
            
            
        end
        
        for num_test=1:n_tests
            
                
            A_id_cell=A_id_cell_RLS_IV_case{num_id_net,num_test};
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
            Constr_list_cell_net{num_case,num_test}=Constr_list;
            
        end
        
       
   
    end
    if not(isempty(Constr_list_cell_net))
        
        num_eval=length(inc_coef2add);
        for idx_edges=1:num_eval
            % matrices used to create an average ranking of each edge over
            % the replicated inference tests for each case
            A_n_rank_1=zeros(n_nodes,n_nodes);
            A_n_rank_2=zeros(n_nodes,n_nodes);
            for num_case=1:size(Constr_list_cell_net,1)
                
                if not(isempty(Constr_list_cell_net(num_case,:)))
                    
                    
                    for idx_id=1:size(Constr_list_cell_net(num_case,:),2)
                        rank_list=Constr_list_cell_net{num_case,idx_id};
                        if not(isempty(rank_list))
                            n_edges=sum(inc_coef2add(1:idx_edges));
                            for idx_n=1:n_edges
                                if num_case==1
                                    A_n_rank_1(rank_list(idx_n,1),rank_list(idx_n,2))=A_n_rank_1(rank_list(idx_n,1),rank_list(idx_n,2))+1;
                                elseif num_case==2
                                    A_n_rank_2(rank_list(idx_n,1),rank_list(idx_n,2))=A_n_rank_2(rank_list(idx_n,1),rank_list(idx_n,2))+1;
                                    
                                end
                            end
                        end
                        
                        
                        
                    end
                end
            end
            
            
            A_n_rank=zeros(n_nodes,n_nodes); % a matrix used to create an average ranking of each edge over the all inference tests

            v_comb=1:num_cases;
            for idx_comb=1:num_cases
                n_comb=nchoosek(v_comb,idx_comb);
                [r_comb, p_comb]=size(n_comb);
                for idx_n_comb=1:r_comb
                    
                    if p_comb==1
                        eval(['A_term=A_n_rank_',num2str(n_comb(idx_n_comb)),';'])
                    elseif p_comb==2
                        
                        eval(['A_term=A_n_rank_',num2str(n_comb(idx_n_comb,1)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,2)),';'])
                        %                             elseif p_comb==3
                        %
                        %                                 eval(['A_term=A_n_rank_',num2str(n_comb(idx_n_comb,1)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,2)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,3)),';'])
                        %                             elseif p_comb==4
                        %
                        %                                 eval(['A_term=A_n_rank_',num2str(n_comb(idx_n_comb,1)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,2)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,3)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,4)),';'])
                        %                             elseif p_comb==5
                        %
                        %                                 eval(['A_term=A_n_rank_',num2str(n_comb(idx_n_comb,1)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,2)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,3)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,4)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,5)),';'])
                        %                             elseif p_comb==6
                        %
                        %                                 eval(['A_term=A_n_rank_',num2str(n_comb(idx_n_comb,1)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,2)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,3)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,4)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,5)),'.*A_n_rank_',num2str(n_comb(idx_n_comb,6)),';'])
                        %
                        
                    end
                    cost_k=nchoosek(num_cases,p_comb);
                    A_term_comb=cost_k* A_term;
                    A_n_rank=A_n_rank+A_term_comb;
                    
                end
                
            end
            A_n_rank_cell_net{num_id_net,idx_edges}=A_n_rank;
            
        end
        
            
    
    
        eval(['load  Data_results_figs_3_4_S1/Data/DREAM4_GoldStandard_InSilico_Size100_',num2str(num_id_net),'.mat'])
        
        
        PPV=zeros(1,num_eval);
        Se=zeros(1,num_eval);
        Sp=zeros(1,num_eval);
        FP=zeros(1,num_eval);
        TN=(n_nodes-1)*n_nodes-nnz(original_network);
        
        
        
        for idx_edg=1:num_eval
            A_n_rank=A_n_rank_cell_net{num_id_net,idx_edg};
            A_id_ones=val_rec_list(A_n_rank,sum(inc_coef2add(1:idx_edg)));
            [Perf_idx, Res]= ppv_sens_spec(original_network,A_id_ones);
            
            PPV(idx_edg)=Perf_idx(2,1);
            Se(idx_edg)=Perf_idx(2,2);
            Sp(idx_edg)=Perf_idx(2,3);
            FP(idx_edg)=size(Res.Dir.False_pos,1);
            
        end
        
        
        
        
        
        aupr_Riva(num_id_net)=eval_aupr(PPV,Se);
        auroc_Riva(num_id_net)=eval_auroc(Se,FP,TN);
    end
                       
end




%%%%%% Plot the performance of RIVA and the other algorithms 

load Data_results_figs_3_4_S1/results_all/Bingo_100
load Data_results_figs_3_4_S1/results_all/dynGenie3_100


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
            name_net='net. 4';
        elseif num_id_net==5
            name_net='net. 5';
    end
%     
    
    f=figure;
    hold on 
    grid on
    plot(auroc_Riva(num_id_net),aupr_Riva(num_id_net),tag_riva,'linewidth',1.5);
       
    plot(Result_dynGenie3100.auroc(num_id_net),Result_dynGenie3100.aupr(num_id_net),tag_genie,'linewidth',1.5);
    plot(Result_Bingo100.auroc(num_id_net),Result_Bingo100.aupr(num_id_net),tag_bingo,'linewidth',1.5);

    xlim([0.5 0.81])
     
    ylim([0 0.25])
    ylabel('AUPR')
    xlabel('AUROC')
    title(name_net)
    set(f,'Position',[10 10 250 250])
    set(gca,'fontsize',12)
    
     if num_id_net==5
        
        legend('RIVA','dynGenie3','Bingo')
    end
    
end
