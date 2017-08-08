% fisher_nonaut_ex2_pert_D_lambda.m written 8-8-17 by JTN to simulate
% nonautonomous Fisher's Equation for example 2 when D(t) = 1 or lambda(t) = 1
% as a means to estimate A or B, where 
% phi(t) = eps1*A*\int D(tau)dtau + eps2*B*\int lambda(tau)dtau

clear all; %clc

%%% perturbing alpha
% par_var = [.4 1 1.5 2];
% pert_var = 'alpha';

%%% perturbing beta
par_var = [3 4 8 12];
pert_var = 'beta';

%%% perturbing gamma
% par_var = [-.15 -.5 -1 -1.5];
% pert_var = 'gamma';

%%% which function is changing over time?
pert_param = 'lambda';

%Construct vectors of independent variables
xn = 301; %number of x points
dt = 1e-3; %time step
t = 0:dt:45;
x = linspace(0,40,xn);
dx = x(2) - x(1);
tn = length(t);

colors = 'brkm';

A_est = zeros(4,1);
B_est = zeros(4,1);

for j = 1:4

    %define activation modulus, signal factor, and response to signal factor
    alpha = 1;
    beta = par_var(j);%4;
    gamma = -1;
    D_large = .25;
    lambda_small = .0125;
    
    [g,sigma,sigma_inv,s,f,int_f_s,psi] = g_sigma_h_example2(alpha,beta,gamma);

    %find x locations where D large
    D_cut = .5;
    lambda_cut = 0.5;

    %these two fixed for speed of 1
    D_small = .25;
    lambda_large = .25;
    
    %nonautonomous diffusion, proliferation definitions

    D_inact = @(t) D_small;
    lambda_inact = @(t) lambda_large;

    D_nonaut = @(t) D_small + (D_large - D_small)*(1-uniform_cdf(0.05,0.35,psi(t)));
    lambda_nonaut = @(t) lambda_large + (lambda_small - lambda_large)*(1-uniform_cdf(0.05,0.35,psi(t)));

    D_pert = @(t) (D_large - D_small)*(1-uniform_cdf(0.05,0.35,psi(t)));
    lambda_pert = @(t) (lambda_small - lambda_large)*(1-uniform_cdf(0.05,0.35,psi(t)));
    
    %initial condition
  
    IC1d = (1-tanh(x-5))/2;   % double(x<5);

    %simulate and time perturbed and unperturbed equations
    tic

        %and nonautonomous
        z_nonaut = RD_sim_nonaut_ex1(D_nonaut,lambda_nonaut,t,x,IC1d);

        z_inact = RD_sim_nonaut_ex1(D_inact,lambda_inact,t,x,IC1d);

    toc

    %plot how far apart the simulations are
 
    figure
    clf
    hold on
    
    LE = cell(4,1);
    LE_inact = cell(4,1);
    LE_u = [0.01 0.1 0.5 0.9];

    %find levels sets for vairous densities
    for i = 1:length(LE_u)
       LE{i} =  leading_edge_calc(z_nonaut(1:100:end,:),x,LE_u(i),0);
       LE_inact{i} = leading_edge_calc(z_inact(1:100:end,:),x,LE_u(i),0);
    end

    %plot the differences in level sets for perturbed, unperturbed
    for i = 1:length(LE_u)
        plot(t(1:100:end),smooth(LE{i}-LE_inact{i}),colors(i))
    end

    

    %numerically integrate D_pert
    lambda_pert_int = zeros(length(t(1:100:end)),1);
    D_pert_int = zeros(length(t(1:100:end)),1);
    
    for i = 1:length(lambda_pert_int)
        lambda_pert_int(i) = 100*dt*trapz(lambda_pert(t(1:100:100*i)));
        D_pert_int(i) = 100*dt*trapz(D_pert(t(1:100:100*i)));
    end

    %plot perturbing function
    if strcmp(pert_param,'lambda')

        plot(t(1:100:end),lambda_pert_int,'--')
        
        est_par = 'B';
        
    elseif strcmp(pert_param,'D')

        plot(t(1:100:end),D_pert_int,'--')
        
        est_par = 'A';
        
    end
        
    xlabel('t')
    ylabel('$\hat{\phi}_{a}(t)$','interpreter','latex')

    if strcmp(pert_param,'lambda')
        legend('a=0.01','a=0.1','a=0.5','a=0.9','\int \lambda(t)dt','location','northeast')
    elseif strcmp(pert_param,'D')
        legend('a=0.01','a=0.1','a=0.5','a=0.9','\int D(t)dt','location','southeast')
    end

    title(['Fisher Homogenization, $\beta$ = ' num2str(par_var(j))] ,'interpreter','latex')

%     exportfig(gcf,['FH_pert_' pert_var '_' num2str(j) '_' est_par '.eps'],'fontsize',2,'color','rgb')

    
    %estimates of A,B
    y=smooth(LE{1}-LE_inact{1})./lambda_pert_int;
    z=smooth(LE{1}-LE_inact{1})./D_pert_int;
    A_est(j) = polyfit(1:200,z(end-199:end)',0);
    B_est(j) = polyfit(1:200,y(end-199:end)',0);
    
end

A_est
B_est