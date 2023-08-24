%% Nonlinear solver for four-species SL model
%% parameters
beta=0.5;
r = 0.84;
mu=0.1;
nu=0.05;
s1=0.01;
s2=0.05;
w0=0.9;
w1=0.4;
t1=0.4;
f1=0.9;
t2=0.4;
f0=0.1;
sigmaF = 0.1;
sigmaR = 0.5;
error_tolerance = 1e-14; % for the nonlinear solver
%% Function definitions
phi = @(f_0, f_1, g, t_2, s_2) f_0 + (f_1-f_0)./(1 + exp(-(g-t_2)/s_2));
fE = @(x, alpha, r, f_0, f_1, t_2, s_2) phi(f_0, f_1, x, t_2, s_2).*(1-x)-alpha.*x.*(1-x).*(1 - (1-x)/r); % ODE nullcline
phi_prime = @(f_0, f_1, x, t_2, s_2) ( (-f_0+f_1)*exp((x+t_2)/s_2) )./(s_2*(exp(x/s_2)+exp(t_2/s_2)).^2);
kernel_FT = @(xi,sigma) exp(-sigma*sigma.*xi.*xi/2);
%% x = [G, S, T], F = 1 - S - G - T
sol_out = [0.2, 0.1, 0.4]; % initial guess
alpha = 0.15:0.01:2;
store_sol = zeros(3,length(alpha));
princ_eig = zeros(2,length(alpha));
for ii = 1:length(alpha)
    f = @(x) [mu*x(2)+nu*x(3)+phi(f0, f1, x(1), t2, s2)*(1-x(1)-x(2)-x(3))-beta*x(1)*x(3)-alpha(ii)*x(1)*(1-x(1)-x(2)-x(3))*(1-(1-x(1)-x(2)-x(3))/r),...
        beta*x(1)*x(3)-phi(w0, w1, x(1), t1, s1)*x(2)-mu*x(2)-alpha(ii)*x(2)*(1-x(1)-x(2)-x(3))*(1-(1-x(1)-x(2)-x(3))/r),...
        phi(w0, w1, x(1), t1, s1)*x(2)-nu*x(3)-alpha(ii)*x(3)*(1-x(1)-x(2)-x(3))*(1-(1-x(1)-x(2)-x(3))/r)];
    options = optimoptions('fsolve','Display','iter','OptimalityTolerance',...
        error_tolerance,'StepTolerance',error_tolerance);
    [sol_out,fval,exitflag,output,jacobian] = fsolve(f,sol_out,options);
    store_sol(1,ii) = sol_out(1);
    store_sol(2,ii) = sol_out(2);
    store_sol(3,ii) = sol_out(3);
    Fcomp = 1 - sol_out(1) - sol_out(2) - sol_out(3);
    % assess stability of FG eigenvalue for spatial instability
    wave_nums = 0:0.01:75;
    [princ_eig(1,ii), princ_eig(2,ii)] = max(alpha(ii)*(1-Fcomp).*(1 - Fcomp/r).*kernel_FT(wave_nums,sigmaF)...
        - (alpha(ii)/r)*Fcomp.*(1-Fcomp).*kernel_FT(wave_nums,sigmaR)...
        - alpha(ii)*Fcomp.*(1-Fcomp/r) - phi(f0, f1, sol_out(1), t2, s2));
end
%% Plotting
plot(alpha,store_sol(1,:),'LineWidth',2);
hold on;
plot(alpha,store_sol(2,:),'LineWidth',2);
plot(alpha,store_sol(3,:),'LineWidth',2);
plot(alpha,1-store_sol(1,:)-store_sol(2,:)-store_sol(3,:),'LineWidth',2);
legend('G','S','T','F');
xline(0.3444,'-.r','LineWidth',2);
xline(0.96535,'-.r','LineWidth',2);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
xlim([min(alpha) max(alpha)]);
ylim([0 1]);