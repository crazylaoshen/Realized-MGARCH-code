function beta = Lassoshoot( y, X, lambda )


 [n, p] = size(X);
 
%  start with ridge regression estimator
 beta = (X'*X + 2*lambda) \ (X'*y);
 
  flag = 0;
    
    % convergence tolerance
    TOL = 1e-6;
    
    while( flag == 0 )
    
        % save current beta
        beta_old = beta;
        
        % optimize elements of beta one by one
        for i = 1:p
            
            % optimize element i of beta
            
            % get ith col of X
            xi = X(:,i);
            
            % get residual excluding ith col
            epsloni = (y - (X*beta - xi*beta(i)));           
            
            % calulate xi'*yi and see where it falls
            Si = (xi'*epsloni); % 1 by 1 scalar
            
            if ( Si < -lambda )
                
                beta(i) = ( Si + lambda )/(xi'*xi);
            
            elseif ( Si > lambda )
            
                beta(i) = ( Si - lambda )/(xi'*xi);
            
            else
                
                beta(i) = 0;
            
            end
            
        end
        
        % check difference between beta and beta_old
        if ( max(abs(beta - beta_old)) <= TOL )
            flag = 1;
        end
                      
    end
 
 
 
 
 
end