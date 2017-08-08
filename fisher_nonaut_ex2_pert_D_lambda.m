% fisher_nonaut_ex2_pert_D_lambda.m written 8-8-17 by JTN to simulate
% nonautonomous Fisher's Equation for example 2 and perform perturbation
% analysis to predict nonautonomous version of Fisher's equation from
% unperturbed equation

clear all; %clc


%Construct vectors of independent variables
xn = 301; %number of x points
dt = 1e-3; %time step
t = 0:dt:45;
x = linspace(0,40,xn);
dx = x(2) - x(1);
tn = length(t);

%previous estimates for A,B
A_est_l = 1;
A_est_u = 1.12;

A_est = 1.06;

B_est_l = .57;
B_est_u = .68;

B_est = 0.6;

% for j = 1:4

    %define activation modulus, signal factor, and response to signal factor
    alpha = 1;
    beta = 4;
    gamma = -.5;
    D_large = .5;
    lambda_small = .0125;
    
    [g,sigma,sigma_inv,s,f,int_f_s,psi] = g_sigma_h_example2(alpha,beta,gamma);

    %find x locations where D large
    D_cut = .5;
    lambda_cut = 0.5;

    %parameter values
    D_small = .25;%D_large*1e-2;
    lambda_large = .25;
    

    %nonautonomous diffusion, proliferation definitions

    D_inact = @(t) D_small;
    lambda_inact = @(t) lambda_large;

    D_nonaut = @(t) D_small + (D_large - D_small)*(1-uniform_cdf(0.05,0.35,psi(t)));
    lambda_nonaut = @(t) lambda_large + (lambda_small - lambda_large)*(1-uniform_cdf(0.05,0.35,psi(t)));

    D_pert = @(t) (D_large - D_small)*(1-uniform_cdf(0.05,0.35,psi(t)));
    lambda_pert = @(t) (lambda_small - lambda_large)*(1-uniform_cdf(0.05,0.35,psi(t)));
    % 
    %initial condition
  
    IC1d = (1-tanh(x-5))/2;   % double(x<5);

    
    %simulate and time perturbed and unperturbed equations
    tic

        z_nonaut = RD_sim_nonaut_ex1(D_nonaut,lambda_nonaut,t,x,IC1d);

        z_inact = RD_sim_nonaut_ex1(D_inact,lambda_inact,t,x,IC1d);

    toc


 %integrate D,lambda

 %numerically integrate D_pert
    lambda_pert_int = zeros(length(t(1:100:end)),1);
    D_pert_int = zeros(length(t(1:100:end)),1);
    
    for i = 1:length(lambda_pert_int)
        lambda_pert_int(i) = 100*dt*trapz(lambda_pert(t(1:100:100*i)));
        D_pert_int(i) = 100*dt*trapz(D_pert(t(1:100:100*i)));
    end

%     phi_r = A_est_u*D_pert_int/(D_large - D_small) + B_est_l*lambda_pert_int/(lambda_small - lambda_large);
%     phi_l = A_est_l*D_pert_int/(D_large - D_small) + B_est_u*lambda_pert_int/(lambda_small - lambda_large);
    
    phi = A_est*D_pert_int*(D_large - D_small) + B_est*lambda_pert_int*(lambda_small - lambda_large);
    
    dt_plot = floor(length(phi)/8);
    
    
    figure
%     clf
    hold on

    for i = 1:8
        plot(x,z_nonaut(i*100*dt_plot,:),'k')

        plot(x-phi(i*dt_plot),z_inact(i*100*dt_plot,:),'r')
        
%         X = [x-phi_l(i*dt_plot) fliplr(x-phi_r(i*dt_plot))];
%         Y = [z_inact(i*100*dt_plot,:) fliplr(z_inact(i*100*dt_plot,:))];
%         
%         h=fill(X,Y,'b','edgecolor','none');
%         set(h,'facealpha',.5)

    end

    plot(t(1:100:end),phi,'b')
    
    legend('Nonaut','Pred')
        
    xlabel('x')
    ylabel('u')

    title(['Fisher Homogenization '] ,'interpreter','latex')

%     exportfig(gcf,['FH_pert_' pert_var '_' num2str(j) '_' est_par '.eps'],'fontsize',2,'color','rgb')
