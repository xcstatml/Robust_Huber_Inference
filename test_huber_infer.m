%huber_inference
clear;
addpath(genpath('./minFunc_2012'));

%% setup parameters and generate data
n=200;
d=5;
sigma=1;
noise_type='g';
noise_para=4;

alpha_list=[0.85,0.9,0.95,0.99];  %0.85:0.01:0.99;
len_alpha=length(alpha_list);
n_boot=1000;

cov_prob_boot_hub=zeros(1, len_alpha);
cov_prob_adap_boot_hub=zeros(1, len_alpha);
n_run=50;
tau_list=zeros(n_run,2);
inter=false; %intercept

tic;
for r=1:n_run
    [X, y, beta_star, noise_sq, noise]=gen_data(n, d, sigma, noise_type, noise_para, inter);
    
    tau_para=1.15;
    weight_type='Gaussian';
    
    %% Huber Bootstrap with initial tau
    
    [beta_hub, tot_loss_hub, loss_hub, tau, beta_OLS, tot_loss_OLS, loss_OLS]...
            =init_huber(X,y,tau_para);
    tau_list(r,1)=tau; % record tau
    [quantile_list_hub, boot_loss_diff_hub]...
        =huber_boot(X,y, n_boot, tau, weight_type, loss_hub, alpha_list);
    loss_diff_hub=calc_huber_loss(X,y,beta_star,tau)-tot_loss_hub;
    cov_prob_boot_hub=cov_prob_boot_hub+(loss_diff_hub<quantile_list_hub);
    
    %% Huber Bootstrap with adaptive tau    
    [adap_beta_hub, adap_tot_loss_hub, adap_loss_hub, adap_tau]=adaptive_tau(X,y,beta_hub,tau);
    tau_list(r,2)=adap_tau; % record tau
    [adap_quantile_list_hub, adap_boot_loss_diff_hub]...
        =huber_boot(X,y, n_boot, adap_tau, weight_type, adap_loss_hub, alpha_list);
    adap_loss_diff_hub=calc_huber_loss(X,y,beta_star,adap_tau)-adap_tot_loss_hub;
    cov_prob_adap_boot_hub=cov_prob_adap_boot_hub+(adap_loss_diff_hub<adap_quantile_list_hub);        
end
toc;

alpha_list    % target level
cov_prob_boot_hub=cov_prob_boot_hub./n_run    % coverage rate for boot_Huber
cov_prob_adap_boot_hub=cov_prob_adap_boot_hub./n_run         % coverage rate for adaptive boot_Huber
