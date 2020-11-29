function [beta_hub,  tot_loss_hub, loss_hub, tau, beta_OLS, tot_loss_OLS, loss_OLS]...
         =init_huber(X,y,tau_para)
% beta_hub: beta parameter from huber
% loss_hub: huber loss for each data
% tot_loss_hub: total huber loss for all data
% tau: parameter tau 
% beta_OLS

    if nargin<3
        tau_para=1;
    end
    
    beta_OLS=X\y;
    res=y-X*beta_OLS;
    loss_OLS=res.^2/2;
    tot_loss_OLS=sum(loss_OLS);
    
    [n,d]=size(X);
    tau=tau_para*(mean(res.^4)*(n/(d+log(n))))^(1/4);
    
    options=struct('METHOD', 'newton', 'display', 'none');
    funObj=@(w)HuberLoss(w,X,y,tau);    
    [beta_hub,tot_loss_hub,exitflag] = minFunc(funObj,zeros(d,1),options);
    if (exitflag<=0)
        warning('Init Huber Regression Fails');
    end      
    
    [~, loss_hub]=calc_huber_loss(X, y, beta_hub, tau);

end