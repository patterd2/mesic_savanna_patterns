%% Function definitions - using Gaussian kernels
GFT = @(x, sigma) exp( -sigma*sigma.*x.*x/2 );
phi = @(f_0, f_1, g, t_2, s_2) f_0 + (f_1-f_0)./(1 + exp(-(g-t_2)/s_2));
fE = @(x, alpha, r, f_0, f_1, t_2, s_2) phi(f_0, f_1, x, t_2, s_2).*(1-x)-alpha.*x.*(1-x).*(1 - (1-x)/r); % ODE nullcline
phi_prime = @(f_0, f_1, x, t_2, s_2) ( (-f_0+f_1)*exp((x+t_2)/s_2) )./(s_2*(exp(x/s_2)+exp(t_2/s_2)).^2);
%% Model parameters

f_0 = 0.1;
f_1 = 0.9;
t_2 = 0.4;
s_2 = 0.05;

sigma_J = 0.1; % seed spread radius
sigma_R = 0.4; % resource competition radius
sigma_W = 0.01; % fire spread radius forest trees

r = 1;
alpha = 2.2;
%r = 0.8;
%alpha = 4.5;

%% Step 0: Find disjoint intervals which contain the equilibria (either one or two equilibria are possible)
dp = 0.00011;
X1 = 0:dp:1; % all the solutions must lie in [0,1]
X2 = X1+dp;
root_locations = fE(X1, alpha, r, f_0, f_1, t_2, s_2).*fE(X2, alpha, r, f_0, f_1, t_2, s_2)<0;
X1 = X1(root_locations);
X2 = X2(root_locations);
num_roots = length(X1);
roots = zeros(num_roots,1);
%% Step 1: Find the value of the homogeneous solution to the nonspatial system
G = zeros(num_roots,1);
for l=1:num_roots
    a = X1(l);
    b = X2(l); % interval which we know contains the zero
    if fE(a,alpha,r, f_0, f_1, t_2, s_2)*fE(b,alpha,r, f_0, f_1, t_2, s_2)>0
        disp('No zero in that interval')
    else
        p = (a + b)/2;
        err = abs( fE(p,alpha,r, f_0, f_1, t_2, s_2) );
        while err > 1e-8
            if fE(a,alpha,r, f_0, f_1, t_2, s_2)*fE(p,alpha,r, f_0, f_1, t_2, s_2)<0
                b = p;
            else
                a = p;
            end
            p = (a + b)/2;
            err = abs( fE(p,alpha,r, f_0, f_1, t_2, s_2) );
        end
    end
    % Asign the equilibrium values
    G(l) = p;
end
%% Step 2: Disregard any unstable equilibria from the ODE - there will sometimes be a saddle between the two stable fixed points within the cusp region
stable = zeros(num_roots,1);
for ii = 1:num_roots
    J_ODE = -(alpha/r)*G(ii)*(1-G(ii)) + alpha*G(ii)*(1-(1-G(ii))/r)-alpha*(1-G(ii))*(1-(1-G(ii))/r)...
        - phi(f_0, f_1, G(ii), t_2, s_2) + phi_prime(f_0, f_1, G(ii), t_2, s_2).*(1-G(ii));
    if J_ODE<0
        stable(ii,1) = 1;
    end
end
%% Step 3: Form the Jacobian at the equilibrium (for each wave number k)
temp_eig = zeros(num_roots,1);
l = 1;
figure(1);
while l < num_roots+1
    if stable(l,1) == 1 % don't check any saddle points, only stable steady-states
        wave_nums = 0:0.025:100;
        J = -(alpha/r)*G(l)*(1-G(l)).*GFT(wave_nums,sigma_R) + alpha*G(l)*(1-(1-G(l))/r).*GFT(wave_nums,sigma_J) -alpha*(1-G(l))*(1-(1-G(l))/r)...
            - phi(f_0, f_1, G(l), t_2, s_2) + phi_prime(f_0, f_1, G(l), t_2, s_2).*(1-G(l)).*GFT(wave_nums,sigma_W);
        plot(wave_nums,J,'.-b');
        hold on;
    end
    l = l+1;
end
xlim([0 30]);
ylim([-1.5 0.5]);
set(gca,'linewidth',2);
set(gca,'FontSize',36);
grid on;
hold on;