function [X, Y, mu, beta_star]=gen_panel_data(n, m, d, sparsity, mu_value,...
                                                  sigma, noise_type, noise_para)
% sparisty is the sparisty level of mu either 0.05 or 0.1
% mu_i = c sqrt(log(m)/n) where c=2,2.5,3,4
    
    if nargin<8
        noise_para=4;
    end 
    if nargin<7
        noise_type='t';
    end
    if nargin<6
        sigma=1;
    end
        
    X=randn(n, d);
    beta_star=-1+2*rand(d, m);
    if noise_type=='t'
        noise=sigma*trnd(noise_para, [n,m]);
    elseif noise_type=='G'
        noise=sigma*randn([n,m]);
    end
    % generate 
    mu=zeros(1,m);
    mu(1:ceil(sparsity*m))=mu_value;
    mu=repmat(mu, n, 1);
    
    Y=mu+X*beta_star+noise;
end