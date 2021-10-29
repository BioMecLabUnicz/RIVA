function f= lambda_GCV(x,X,H,X_tr,y)
% define the function f based on generalized cross validation for obtaining 
% a good estimate of x (lambda, i.e. alpha in the manuscript,  ridge parameter) 


n_tps=size(X,1);

I=eye(size(H,1));

A_lambda=X/(H+n_tps*x*I)*X_tr;
%I_a=sparse(eye(size(A_lambda,1)));
f=1/n_tps*norm((sparse(eye(size(A_lambda,1)))-A_lambda)*y)^2/(1/n_tps*trace(sparse(eye(size(A_lambda,1)))-A_lambda))^2;
