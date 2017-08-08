%watch video of unpert v pert sims

%Construct vectors of independent variables
mn = 81; %number of m points
xn = 601; %number of x points
total = mn*xn;
dt = 1e-3; %time step
t = 0:dt:45;
m = linspace(0,1,mn);
dm = m(2) - m(1);
x = linspace(0,40,xn);
dx = x(2) - x(1);
[X,M] = meshgrid(x,m);
tn = length(t);
m_fine = [linspace(0,0.1,100) linspace(0.1,0.9,100) linspace(0.9,1,100)];

%define activation modulus, signal factor, and response to signal factor
alpha = 1;
% beta = 2.005;
beta = 4;
gamma = -.15;
D_large = .25;
lambda_small = .0125;


[g,sigma,sigma_inv,s,f,int_f_s,psi] = g_sigma_h_example2(alpha,beta,gamma);

IC_1_d_m = IC_uniform(.05,.35);

Soln = @(t,s) g(sigma_inv(-int_f_s(t),s))./(g(s)).*IC_1_d_m(sigma_inv(-int_f_s(t),s));


%find x locations where D large
D_cut = .5;
lambda_cut = 0.5;

%parameter values

D_small = .25;%D_large*1e-2;

lambda_large = .25;
%0.3;%lambda_large*1e-2;

%nonautonomous diffusion, proliferation

D_inact = @(t) D_small;
lambda_inact = @(t) lambda_large;

D_nonaut = @(t) D_small + (D_large - D_small)*(1-uniform_cdf(0.05,0.35,psi(t)));
lambda_nonaut = @(t) lambda_large + (lambda_small - lambda_large)*(1-uniform_cdf(0.05,0.35,psi(t)));

D_pert = @(t) (D_large - D_small)*(1-uniform_cdf(0.05,0.35,psi(t)));
lambda_pert = @(t) (lambda_small - lambda_large)*(1-uniform_cdf(0.05,0.35,psi(t)));
% 
%initial condition

IC1d = (1-tanh(x-5))/2;   % double(x<5);

tic

%and nonautonomous
z_nonaut = RD_sim_nonaut_ex1(D_nonaut,lambda_nonaut,t,x,IC1d);

z_inact = RD_sim_nonaut_ex1(D_inact,lambda_inact,t,x,IC1d);

toc

figure

for i = 1:1000:tn
    
    hold off
    plot(x,z_nonaut(i,:))
    hold on
    plot(x,z_inact(i,:))
    
    title(['t = ' num2str(t(i))])
    xlabel('x')
    ylabel('u')
    
    axis([0 40 0 1])
    
    pause(.125)
end