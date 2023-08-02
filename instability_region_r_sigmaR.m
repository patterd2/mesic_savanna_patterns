%% Function definitions
phi = @(f_0, f_1, g, t_2, s_2) f_0 + (f_1-f_0)./(1 + exp(-(g-t_2)/s_2));
fE = @(x, alpha, r, f_0, f_1, t_2, s_2) phi(f_0, f_1, x, t_2, s_2).*(1-x)-alpha.*x.*(1-x).*(1 - (1-x)/r); % ODE nullcline
phi_prime = @(f_0, f_1, x, t_2, s_2) ( (-f_0+f_1)*exp((x+t_2)/s_2) )./(s_2*(exp(x/s_2)+exp(t_2/s_2)).^2);
kernel_FT = @(xi,sigma) exp(-sigma*sigma.*xi.*xi/2);
%% Model parameters

f_0 = 0.1;
f_1 = 0.9;
t_2 = 0.4;
s_2 = 0.05;

sigmaW = 0.01;
sigmaF = 0.1;

sigmaR = 0.005:0.001:0.15;
r = 0.6:0.001:1;
alpha = 3.75;
conditions = -10*ones(length(r),length(sigmaR));
tic;
for tt = 1:length(sigmaR)
    for kk = 1:length(r)
        %% Step 0: Find disjoint intervals which contain the equilibria (either one or two equilibria are possible)
        dp = 0.00051;
        X1 = -0.01:dp:1.01; % all the solutions must lie in [0,1]
        X2 = X1+dp;
        root_locations = fE(X1, alpha, r(kk), f_0, f_1, t_2, s_2).*fE(X2, alpha, r(kk), f_0, f_1, t_2, s_2)<0;
        X1 = X1(root_locations);
        X2 = X2(root_locations);
        num_roots = length(X1);
        roots = zeros(num_roots,1);
        %% Step 1: Find the value of the homogeneous solution to the nonspatial system
        G = zeros(num_roots,1);
        for l=1:num_roots
            a = X1(l);
            b = X2(l); % interval which we know contains the zero
            if fE(a,alpha,r(kk), f_0, f_1, t_2, s_2)*fE(b,alpha,r(kk), f_0, f_1, t_2, s_2)>0
                disp('No zero in that interval')
            else
                p = (a + b)/2;
                err = abs( fE(p,alpha,r(kk), f_0, f_1, t_2, s_2) );
                while err > 1e-10
                    if fE(a,alpha,r(kk), f_0, f_1, t_2, s_2)*fE(p,alpha,r(kk), f_0, f_1, t_2, s_2)<0
                        b = p;
                    else
                        a = p;
                    end
                    p = (a + b)/2;
                    err = abs( fE(p,alpha,r(kk), f_0, f_1, t_2, s_2) );
                end
            end
            % Asign the equilibrium values
            G(l) = p;
        end
        %% Step 2: Disregard any unstable equilibria from the ODE - there will sometimes be a saddle between the two stable fixed points within the cusp region
        stable = zeros(num_roots,1);
        for ii = 1:num_roots
            J_ODE = -(alpha/r(kk))*G(ii)*(1-G(ii)) + alpha*G(ii)*(1-(1-G(ii))/r(kk))-alpha*(1-G(ii))*(1-(1-G(ii))/r(kk))...
                - phi(f_0, f_1, G(ii), t_2, s_2) + phi_prime(f_0, f_1, G(ii), t_2, s_2).*(1-G(ii));
            if J_ODE <= 0
                stable(ii,1) = 1;
            end
        end
        G = G(stable==1); % throw away unstable fixed points
        F = 1 - G;
        %% Step 3: Compute the sufficient condition for pattern formation
        wave_nums = 0:0.01:25;
        if isempty(G)==0
            for index = 1:length(G)
                conditions(kk,tt) = max(conditions(kk,tt), max(phi_prime(f_0, f_1, 1-F(index), t_2, s_2).*F(index).*kernel_FT(wave_nums,sigmaW) +...
                    - (alpha*F(index).*(1-F(index))/r(kk)).*kernel_FT(wave_nums,sigmaR(tt))...
                    + alpha*(1-F(index)).*(1-F(index)/r(kk)).*kernel_FT(wave_nums,sigmaF) ...
                    - alpha*F(index).*(1-F(index)/r(kk))...
                    - phi(f_0, f_1, 1-F(index), t_2, s_2)) );
            end
        end
        if isempty(G)
            conditions(kk,tt) = 0;
        end
    end
end
%% plotting 
imagesc(r,sigmaR,conditions');
colormap(bluewhitered);
ylabel('\sigma_R');
xlabel('r');
set(gca,'YDir','normal');
set(gca,'linewidth',2);
set(gca,'FontSize',36);
xlim([min(r) max(r)]);
ylim([min(sigmaR) max(sigmaR)]);
xticks([0.6 0.7 0.8 0.9 1]);
yticks([0 0.05 0.1 0.15]);

%%
toc;

