function [dec_seq, FDP, power]=BH(p_seq, sig_seq, alpha)

    N=length(sig_seq);
    [p_sorted, sort_ids]=sort(p_seq);
    thresh=(1:N)*alpha/N;
    rej=p_sorted<=thresh;
    max_id=find(rej,1,'last'); 
    dec_seq=zeros(1,N);
    dec_seq(sort_ids(1:max_id))=1;
    
    false_rej=sum(dec_seq.*(1-sig_seq));
    rej=sum(dec_seq);
    FDP=false_rej/max(sum(dec_seq),1);
    power=(rej-false_rej)/max(sum(sig_seq),1);
     
   
end