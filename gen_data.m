function [X, y, beta_star, noise_sq, noise]=gen_data(n, d, sigma, noise_type, noise_para, inter)
    
    if nargin<6
        inter=false;
    end
    
    if nargin<5
        noise_para=2;
    end 
    if nargin<4
        noise_type='t';
    end
    if nargin<3
        sigma=1;
    end
     
    if inter % including an intercept
        X=[ones(n,1), randn(n,d-1)];
        beta_star=ones(d,1);
        beta_star(2:d)=rand(d-1,1)<0.5;
        beta_star(2:d)=(beta_star(2:d)-0.5)*2;
    else
        X=randn(n,d);
        beta_star=rand(d,1)<0.5;
        beta_star=1.5*(beta_star-0.5)*2;
    end
    
    %X=randn(n, d-1);
    %X=[ones(n,1),X];    
    %beta_star=linspace(0, 1, d)';
    %beta_star=rand(d-1,1)<0.5;
    %beta_star=(beta_star-0.5)*2;
    %beta_star=[2;beta_star];
    
    if noise_type=='t'
        noise=sigma*trnd(noise_para, [n,1]);
    elseif noise_type=='g'
        noise=noise_para*randn([n,1]);
    elseif  strcmp(noise_type,'log') %lognormal(mu,sig)
        mu=0;
        sig=noise_para;
        std_noise=1;%sqrt((exp(sig^2)-1)*exp(sig^2));
        noise=sigma*(lognrnd(mu,sig,[n,1])-exp(sig^2/2)*ones(n,1))/std_noise;        
    elseif strcmp(noise_type, 'par')
        mean_noise=noise_para/(noise_para-1);
        noise=gprnd(1/noise_para,1/noise_para,1,n,1);
        std_noise=1; %sqrt(noise_para^3/(noise_para-1)^2/(noise_para-2));
        %noise=sqrt(2)*(sigma*(wt*randn([n,1])+(1-wt)*(noise-mean_noise*ones(n,1))/std_noise));
        noise=sigma*(noise-mean_noise*ones(n,1))/std_noise;
    end
    noise_sq=norm(noise,2)^2/2;
    y=X*beta_star+noise;   
end