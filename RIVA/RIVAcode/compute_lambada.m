
function lambda_min=compute_lambada(min_eig_H,max_eig_H)
% compute a ridge regression parameter

if min_eig_H<max_eig_H && min_eig_H>=max_eig_H/3
    lambda_min=max_eig_H/2;
elseif min_eig_H<max_eig_H/3 && min_eig_H>=max_eig_H/5
    lambda_min=max_eig_H/4;
elseif min_eig_H<max_eig_H/5 && min_eig_H>=max_eig_H/7
    lambda_min=max_eig_H/6;
elseif min_eig_H<max_eig_H/7 && min_eig_H>=max_eig_H/9
    lambda_min=max_eig_H/8;
elseif min_eig_H<max_eig_H/9 && min_eig_H>=max_eig_H/11
    lambda_min=max_eig_H/10;
elseif min_eig_H<max_eig_H/11 && min_eig_H>=max_eig_H/20
    lambda_min=max_eig_H/11;
elseif min_eig_H<max_eig_H/20 && min_eig_H>=max_eig_H/30
    lambda_min=max_eig_H/20;
elseif min_eig_H<max_eig_H/30 && min_eig_H>=max_eig_H/40
    lambda_min=max_eig_H/30;
elseif min_eig_H<max_eig_H/40 && min_eig_H>=max_eig_H/50
    lambda_min=max_eig_H/40;    
elseif min_eig_H<max_eig_H/50 && min_eig_H>=max_eig_H/60
    lambda_min=max_eig_H/50;
elseif min_eig_H<max_eig_H/60 && min_eig_H>=max_eig_H/70
    lambda_min=max_eig_H/60;
elseif min_eig_H<max_eig_H/70 && min_eig_H>=max_eig_H/80
    lambda_min=max_eig_H/70;
elseif min_eig_H<max_eig_H/80 && min_eig_H>=max_eig_H/90
    lambda_min=max_eig_H/80;
elseif min_eig_H<max_eig_H/90 && min_eig_H>=max_eig_H/100
    lambda_min=max_eig_H/90;
else
    lambda_min=max_eig_H/100;
end

% if lambda_min>10
%     
%     lambda_min=max_eig_H/min_eig_H;
% end
% 
if lambda_min>100
    lambda_min=lambda_min/10;
end