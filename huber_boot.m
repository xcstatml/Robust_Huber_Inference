function [quantile_list, boot_loss_diff, boot_value, org_value, exitflag]...
          =huber_boot(X, y, n_boot, tau, weight_type, loss_hub, alpha_list)

% Data X in n by d, Y in n by 1.
% n_boot: number of bootstrap samples
% tau: Huber threshold
% weight_type: standardized i.i.d. random weights with mean and variance 1
% loss_hub: huber loss for each data
% alpha_list: confidence level 
      
    [n ,d] =size(X);
    boot_loss_diff=zeros(n_boot, 1);
    exitflag=zeros(n_boot, 1);
    boot_value=zeros(n_boot, 1);
    org_value=zeros(n_boot, 1);
    
    
    for boot_idx=1:n_boot
        if strcmpi(weight_type, 'Gaussian')
            weight=randn(n,1)+1;
            %options=struct('Method', 'lbfgs', 'display', 'none');            
        elseif strcmpi(weight_type, 'exp')
            weight=exprnd(1, [n,1]);            
        elseif strcmpi(weight_type, 'Ber')
            weight=2*binornd(1,0.5,[n,1]);                        
        end
        options=struct('METHOD', 'newton', 'display', 'none');
        funObj=@(w)WeightHuberLoss(w,X,y,tau,weight);
        [~,boot_value(boot_idx),exitflag(boot_idx)] = minFunc(funObj,zeros(d,1),options);
        if exitflag(boot_idx)<0
            warning(['Error in Bootstrap Iter ', num2str(boot_idx)]);
        end
        
        org_value(boot_idx)=dot(weight, loss_hub);
        boot_loss_diff(boot_idx) = org_value(boot_idx)-boot_value(boot_idx);        
    end
   
    quantile_list=quantile(boot_loss_diff, alpha_list);
end