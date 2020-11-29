function [mu_hub,tau,beta_hub]=init_panel_huber(X,Y,tau_para)
% mu_hub: mu parameter from huber
% beta_hub: beta parameter from huber
% tau: parameter tau 

    if nargin<3
        tau_para=1.2;
    end
    
    [d]=size(X,2);
    [n,m]=size(Y);
    
    
    options=struct('Method', 'netwon', 'display', 'none');
    beta_hub=zeros(d,m);   
    tau=zeros(m,1);
    
    for k=1:m
        %% compute tau parameter
        beta_OLS=X\Y(:,k);
        res=Y(:,k)-X*beta_OLS;
        tau(k)=tau_para*(mean(res.^4)*(n/(d+log(n))))^(1/4);
        
        %% compute huber regression
        funObj=@(w)HuberLoss(w,X,Y(:,k),tau(k));    
        [beta_hub(:,k),tot_hub_loss,exitflag] = minFunc(funObj,zeros(d,1),options);
        if (exitflag<=0)
            warning('Init Huber Regression Fails');
        end      
    end
    mu_hub=beta_hub(1,:);
end