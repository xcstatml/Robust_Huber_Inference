function [tot_loss, loss]=calc_huber_loss(X,y,beta,tau)

    n=length(y);
    res=y-X*beta;
    closeInd = abs(res) <= tau;
    loss=zeros(n,1);
    loss(closeInd)=(1/2)*res(closeInd).^2;
    loss(~closeInd)=tau*abs(res(~closeInd)) - (1/2)*tau^2;    
    tot_loss=sum(loss);
    
end