function [beta_hub, tot_loss_hub, loss_hub, tau]=adaptive_tau(X,y, beta_hub, tau, thresh)
    
    if nargin<6
        thresh=1e-5;
    end
    
    [n,d]=size(X);
    options=struct('METHOD', 'newton', 'display', 'none');
    iter=1;
    maxiter=100;
    while iter<maxiter
        %% update for the tau
        res_hub=y-X*beta_hub;
        res4=res_hub.^4;
        tau_search=@(x) fixed_equ_value(x, res4, n, d);
        [tau4, fval, flag]=fzero(tau_search, [min(res4), sum(res4)]);
        if flag<0 %fail case
            tau4=sum(res4)/(n-d);
        end        
        tau_new =tau4^(1/4);
        
        
        %% update for beta         
         funObj=@(w)HuberLoss(w,X,y,tau_new);    
        [beta_new] = minFunc(funObj,zeros(d,1),options);
        
        if max(norm(beta_new-beta_hub,2),abs(tau_new-tau))<thresh
            tau=tau_new;
            beta_hub=beta_new;
            break;
        else
            tau=tau_new;
            beta_hub=beta_new;
            iter=iter+1;            
        end
    end
    
     [tot_loss_hub, loss_hub]=calc_huber_loss(X, y, beta_hub, tau);
    
end

function f=fixed_equ_value(x, res4, n, d)
    f=sum(min(res4,x))/(x*(n-d)) - (log(n)+d)/n;
end