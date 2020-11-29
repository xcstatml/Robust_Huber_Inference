function [p_value, mu_boot]...
          =huber_panel_boot(X, Y, n_boot, tau, weight_type, mu_hub)

% Data X in n by d, Y in n by 1.
% n_boot: number of bootstrap samples
% tau: Huber threshold
% weight_type: standardized i.i.d. random weights with mean and variance 1
% mu_hub: mu parameter from huber regression
      
    [n ,d] =size(X);
    m=size(Y,2);
            
    mu_boot=zeros(m, n_boot);
    p_value=zeros(1,m);
    
    options=struct('Method', 'newton', 'display', 'none');
    
    for boot_idx=1:n_boot
        if strcmpi(weight_type, 'Gaussian')
            weight=randn(n,m)+1;           
        elseif strcmpi(weight_type, 'exp')
            weight=exprnd(1, [n,m]);            
        elseif strcmpi(weight_type, 'Ber')
            weight=2*binornd(1,0.5,[n,m]);                        
        end
        for k=1:m
            funObj=@(w)WeightHuberLoss(w,X,Y(:,k),tau(k),weight(:,k));
            [beta_boot, boot_value, exitflag] = minFunc(funObj,zeros(d,1),options);
            if exitflag<0
                warning(['Error in Bootstrap Iter ', num2str(boot_idx), ' , ', num2str(m)]);
            end
            mu_boot(k, boot_idx)=beta_boot(1);            
            p_value(k)=p_value(k)+(abs(mu_boot(k, boot_idx)-mu_hub(k))>=abs(mu_hub(k)));
        end               
    end
    p_value=p_value/(n_boot+1);
 
end