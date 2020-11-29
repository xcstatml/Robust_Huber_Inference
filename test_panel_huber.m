%huber_inference
clear;
addpath(genpath('./minFunc_2012'));

%% setup parameters and generate data
n=100;
d=3;
m=1000;
sigma=1;
noise_type='t';
noise_para=3.5;

alpha_list=[0.05, 0.1, 0.15, 0.2, 0.25];
len_alpha=length(alpha_list);

mu_value=3*sqrt(2*log(m)/n);
sparsity=0.05;

tau_para=1.15;

n_boot=500;
n_run=50;

%% Multiple Huber Inference
FDP=zeros(n_run, len_alpha);
power=zeros(n_run, len_alpha);
dec_seq=cell(n_run, len_alpha);

for r=1:n_run
    tic;
    [X, Y, mu, beta_star]=gen_panel_data(n, m, d, sparsity, mu_value,...
                                         sigma, noise_type, noise_para);  
    sig_seq=(mu(1,:)>eps);
    
    aug_X=[ones(n,1), X];
    %% initial estimate
    [mu_hub, tau, beta_hub]=init_panel_huber(aug_X,Y,tau_para);
    weight_type='Gaussian';
    
    [p_value, mu_boot]...
          =huber_panel_boot(aug_X, Y, n_boot, tau, weight_type, mu_hub);  
    for alpha_idx=1:len_alpha
        alpha=alpha_list(alpha_idx);
        [dec_seq{r, alpha_idx}, FDP(r,alpha_idx), power(r, alpha_idx)]=...
                                            BH(p_value, sig_seq, alpha);
    end
    time=toc;
    fprintf('r=%d, time=%.2f\n', r, time);
end

alpha_list
avg_FDP=mean(FDP, 1)
avg_power=mean(power, 1)
