    
function [Perf_idx,Res_Struct] = ppv_sens_spec(Orig_mtx,Rec_mtx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a network of n nodes, whose dynamics are described by the evolution
% of a linear time-invariant system dx(t)/dt = A x(t) + B u(t), with n 
% states, this function computes the following reconstruction performance 
% indexes:
%
% PPV (positive predictive value), sensitivity and specificity.
%
%                                 number of true positives  
% PPV           =    ----------------------------------------------------
%                    number of true positives + number of false positives
%
%                                 number of true positives  
% sensitivity   =    ----------------------------------------------------
%                    number of true positives + number of false negatives
%
%                                 number of true negatives  
% specificity   =    ----------------------------------------------------
%                    number of true negatives + number of false positives
%
%
% INPUT VARIABLES
% 
% Orig_mtx:     Original connectivity matrix / network
% 
% Rec_mtx:      Reconstruced connectivity matrix / network
%
% OUTPUT VARIABLES
%
% Perf_idx: has the following structure
%
%       | PPV | SENSITIVITY | SPECIFICITY |
% ----------------------------------------|
% UNDIR |  *  |      *      |     *       |
% ----------------------------------------|
% DIR   |  *  |      *      |     *       |
% ----------------------------------------|
% SIGN  |  *  |      *      |     *       |
% -----------------------------------------
%
%                      first column -->  ppv value                                     
%                      second column --> sensitivity value 
%                      third column --> specificity value 
%                      first row --> performance evaluation without taking
%                                    into account the direction of the arcs                                    
%                                    of network / matrix connectivity
%                                    (UNDIRECTED NETWORK)
%                      second row --> performance evaluation taking
%                                     into account the direction of the
%                                     arcs of network / matrix connectivity
%                                     (DIRECTED NETWORK)
%                      third row --> performance evaluation taking
%                                    into account the direction and the 
%                                    influence type (promoting or 
%                                    inhibiting)  of the arcs of
%                                    network / matrix connectivity
%                                    (SIGNED NETWORK)
% Res_Struct is a structure containing the true positives, false positives,
% true negatives and false negatives for each case (UNDIRECTED, DIRECTED
% and SIGNED network)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n_row,n_col]=size(Rec_mtx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORMANCE EVALUATION : UNDIRECTED NETWORK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% true positives counter
true_positive_undir=0;
% true negatives counter
true_negative_undir=0;
% false positives counter
false_positive_undir=0;
% false negatives counter
false_negative_undir=0;

% True Positives Location Matrix (Undirected network)
Res_Struct.Und.True_pos=[];
% True Negatives Location Matrix
Res_Struct.Und.True_neg=[];
% False Positives Location Matrix
Res_Struct.Und.False_pos=[];
% False Negatives Location Matrix
Res_Struct.Und.False_neg=[];


for idx1=1:n_row
    for idx2=idx1+1:n_col
        % true positives occurrence
        if ((Orig_mtx(idx1,idx2) || Orig_mtx(idx2,idx1)) ...
                && (Rec_mtx(idx1,idx2) || Rec_mtx(idx2,idx1)))
            true_positive_undir=true_positive_undir+1;
            Res_Struct.Und.True_pos=[Res_Struct.Und.True_pos ; [idx1,idx2]];
        end
        % true negatives occurrence
        if ((not(Orig_mtx(idx1,idx2)) && not(Orig_mtx(idx2,idx1))) ...
                && (not(Rec_mtx(idx1,idx2)) && not(Rec_mtx(idx2,idx1))))
            true_negative_undir=true_negative_undir+1;
            Res_Struct.Und.True_neg=[Res_Struct.Und.True_neg; [idx1,idx2]];
        end
        % false positives occurrence
        if ((not(Orig_mtx(idx1,idx2)) && not(Orig_mtx(idx2,idx1))) ...
                && (Rec_mtx(idx1,idx2) || Rec_mtx(idx2,idx1)))
            false_positive_undir=false_positive_undir+1;
            Res_Struct.Und.False_pos=[Res_Struct.Und.False_pos; [idx1,idx2]];
        end
        % false negatives occurrence
        if ((Orig_mtx(idx1,idx2) || Orig_mtx(idx2,idx1)) ...
                && (not(Rec_mtx(idx1,idx2)) && not(Rec_mtx(idx2,idx1))))
            false_negative_undir=false_negative_undir+1;
            Res_Struct.Und.False_neg=[Res_Struct.Und.False_neg; [idx1,idx2]];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ppv, sensitivity and specificity evaluation without assuming the 
% direction of the arcs of the network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ppv evalutation
PPV_undir=true_positive_undir/(true_positive_undir+false_positive_undir);
% sensitivity evaluation
sensitivity_undir=true_positive_undir/(true_positive_undir+false_negative_undir);
% specificity evaluation
specificity_undir=true_negative_undir/(false_positive_undir+true_negative_undir);
% row vector containing ppv, sensitivity and specificity values
perf_ind_undir=[PPV_undir,sensitivity_undir,specificity_undir];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORMANCE EVALUATION : DIRECTED-SIGNED NETWORK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% true positives counter
true_positive_sign=0;
% true negatives counter
true_negative=0;
% false positives counter
false_positive_dir=0;
% false negatives counter
false_negative_dir=0;
% sign mismatches counter
sign_mismatch=0;

% True Positives Location Lists (Directed network)
Res_Struct.Dir.True_pos=[];
% True Negatives Location List
Res_Struct.Dir.True_neg=[];
% False Positives Location List
Res_Struct.Dir.False_pos=[];
% False Negatives Location List
Res_Struct.Dir.False_neg=[];

% True Positives Location Lists (Signed-Directed network)
Res_Struct.Sign.True_pos=[];
% True Negatives Location List
Res_Struct.Sign.True_neg=[];
% False Positives Location List
Res_Struct.Sign.False_pos=[];
% False Negatives Location List
Res_Struct.Sign.False_neg=[];
% Sign Mismatches Location List
Res_Struct.Sign.Mism=[];

for idx1=1:n_row
    for idx2=[1:idx1-1,idx1+1:n_col]
        if(Rec_mtx(idx1,idx2) && Orig_mtx(idx1,idx2))
            % true positive and sign mismatch occurrence assuming the
            % influence type of the arcs
            if ((Rec_mtx(idx1,idx2)*Orig_mtx(idx1,idx2))>0) % same signs
                true_positive_sign=true_positive_sign+1;
                Res_Struct.Sign.True_pos=[Res_Struct.Sign.True_pos; [idx1,idx2]];                
            else % different signs
                sign_mismatch=sign_mismatch+1;
                Res_Struct.Sign.Mism=[Res_Struct.Sign.Mism; [idx1,idx2]];              
            end
            % false positive occurrence
        elseif (Rec_mtx(idx1,idx2) && not(Orig_mtx(idx1,idx2)))
            false_positive_dir=false_positive_dir+1;
            Res_Struct.Dir.False_pos=[Res_Struct.Dir.False_pos; [idx1,idx2]];
            % false negative occurrence
        elseif (not(Rec_mtx(idx1,idx2)) && Orig_mtx(idx1,idx2))
            false_negative_dir=false_negative_dir+1;
            Res_Struct.Dir.False_neg=[Res_Struct.Dir.False_neg; [idx1,idx2]];             
        else
            % true negative occurrence
            true_negative=true_negative+1;
            Res_Struct.Dir.True_neg=[Res_Struct.Dir.True_neg; [idx1,idx2]];            
        end
    end
end

% False negatives occurrence assuming the influence type of arcs
% (signed network)
false_negative_sign = false_negative_dir+sign_mismatch;
Res_Struct.Sign.False_neg=[Res_Struct.Dir.False_neg;Res_Struct.Sign.Mism];

% True negatives are the same for the DIR and SIGN cases
Res_Struct.Sign.True_neg=Res_Struct.Dir.True_neg;

% True positives occurrence without assuming the influence type of arcs
% (directed network)
true_positive_dir = true_positive_sign+sign_mismatch;
Res_Struct.Dir.True_pos=[Res_Struct.Sign.True_pos; Res_Struct.Sign.Mism];

% False positives occurrence assuming the influence type of arcs
% (signed network)
false_positive_sign = false_positive_dir+sign_mismatch;
Res_Struct.Sign.False_pos=[Res_Struct.Dir.False_pos; Res_Struct.Sign.Mism];            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ppv, sensitivity and specificity evaluation assuming only the direction
% of the arcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ppv evaluation
PPV_dir=true_positive_dir/(true_positive_dir+false_positive_dir);
% sensitivity evaluation
sensitivity_dir=true_positive_dir/(true_positive_dir+false_negative_dir);
% specificity evaluation
specificity_dir=true_negative/(false_positive_dir+true_negative);
% row vector containing ppv, sensitivity and specificity values
perf_ind_dir=[PPV_dir,sensitivity_dir,specificity_dir];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ppv, sensitivity and specificity evaluation assuming the direction and
% the influence type of the arcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ppv evaluation
PPV_sign=true_positive_sign/(true_positive_sign+false_positive_sign);
% sensitivity evaluation
sensitivity_sign=true_positive_sign/(true_positive_sign+false_negative_sign);
% specificity evaluation
specificity_sign=true_negative/(false_positive_sign+true_negative);
% row vector containing ppv, sensitivity and specificity values
perf_ind_sign=[PPV_sign,sensitivity_sign,specificity_sign];

% matrix containing performance indexes
Perf_idx =[perf_ind_undir;perf_ind_dir;perf_ind_sign];