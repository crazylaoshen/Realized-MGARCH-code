function lambda = LassoLambda( y, X, kfold, lambda_vec )
%   Estimate L1 regularization parameter in LASSO via cross validation.
%   Algorithm for solving the Lasso problem:
%           0.5 * (y - X*beta)'*(y - X*beta) + lambda * ||beta||_1                                       
%   where ||beta||_1 is the L_1 norm i.e., ||beta||_1 = sum(abs( beta ))



 % get size of X
    [n, p] = size(X);


%   partition input data for cross-validation
    
    CVO = cvpartition( n, 'kfold' , kfold );


MSEerror = zeros( length(lambda_vec), 1 );

    for k = 1:length( lambda_vec )
        
      % vector to store MSE error for each fold  
      err = zeros(CVO.NumTestSets,1);
      
      for i = 1:CVO.NumTestSets
          
          % get ith training set
          trIdx = CVO.training(i);
          
          % get ith test set
          teIdx = CVO.test(i);
          
          % train LASSO using training set
          beta = Lassoshoot( y(trIdx), X(trIdx,:), lambda_vec(k) );

          % test using testing set          
          ypred = X(teIdx,:) * beta;
                    
          % calculate MSE error            
          temp = y(teIdx) - ypred;
                    
          % clear ypred
          clear ypred;
          
          % calculate average error for the ith test set
          err(i) = (temp'*temp)/length(temp);  
                    
      end
      
      % calculate the mean error over all test sets
      MSEerror(k) = mean(err);
    end
    
        % find the smallest MSE
    lambda = lambda_vec( MSEerror == min(MSEerror) );
    
    


end