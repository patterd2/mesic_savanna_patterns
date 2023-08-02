%% Function definitions
phi = @(f_0, f_1, g, t_2, s_2) f_0 + (f_1-f_0)./(1 + exp(-(g-t_2)/s_2));
fE = @(x, alpha, r, f_0, f_1, t_2, s_2) phi(f_0, f_1, x, t_2, s_2).*(1-x)-alpha.*x.*(1-x).*(1 - (1-x)/r); % ODE nullcline
phi_prime = @(f_0, f_1, x, t_2, s_2) ( (-f_0+f_1)*exp((x+t_2)/s_2) )./(s_2*(exp(x/s_2)+exp(t_2/s_2)).^2);
%% Model parameters

f_0 = 0.1;
f_1 = 0.9;
t_2 = 0.4;
s_2 = 0.05;

r = 0.01:0.001:1;

alpha = 0.901:0.001:6;
conditions = zeros(length(alpha),length(r));
tic;
for jj = 1:length(alpha)
    for kk = 1:length(r)
        %% Step 0: Find disjoint intervals which contain the equilibria (either one or two equilibria are possible)
        dp = 0.0011;
        X1 = -0.01:dp:1.01; % all the solutions must lie in [0,1]
        X2 = X1+dp;
        root_locations = fE(X1, alpha(jj), r(kk), f_0, f_1, t_2, s_2).*fE(X2, alpha(jj), r(kk), f_0, f_1, t_2, s_2)<0;
        X1 = X1(root_locations);
        X2 = X2(root_locations);
        num_roots = length(X1);
        roots = zeros(num_roots,1);
        %% Step 1: Find the value of the homogeneous solution to the nonspatial system
        G = zeros(num_roots,1);
        for l=1:num_roots
            a = X1(l);
            b = X2(l); % interval which we know contains the zero
            if fE(a,alpha(jj),r(kk), f_0, f_1, t_2, s_2)*fE(b,alpha(jj),r(kk), f_0, f_1, t_2, s_2)>0
                disp('No zero in that interval')
            else
                p = (a + b)/2;
                err = abs( fE(p,alpha(jj),r(kk), f_0, f_1, t_2, s_2) );
                while err > 1e-8
                    if fE(a,alpha(jj),r(kk), f_0, f_1, t_2, s_2)*fE(p,alpha(jj),r(kk), f_0, f_1, t_2, s_2)<0
                        b = p;
                    else
                        a = p;
                    end
                    p = (a + b)/2;
                    err = abs( fE(p,alpha(jj),r(kk), f_0, f_1, t_2, s_2) );
                end
            end
            % Asign the equilibrium values
            G(l) = p;
        end
        %% Step 2: Disregard any unstable equilibria from the ODE - there will sometimes be a saddle between the two stable fixed points within the cusp region
        stable = zeros(num_roots,1);
        for ii = 1:num_roots
            J_ODE = -(alpha(jj)/r(kk))*G(ii)*(1-G(ii)) + alpha(jj)*G(ii)*(1-(1-G(ii))/r(kk))-alpha(jj)*(1-G(ii))*(1-(1-G(ii))/r(kk))...
                - phi(f_0, f_1, G(ii), t_2, s_2) + phi_prime(f_0, f_1, G(ii), t_2, s_2).*(1-G(ii));
            if J_ODE<0
                stable(ii,1) = 1;
            end
        end
        G = G(stable==1); % throw away unstable fixed points
        %% Step 3: Compute the sufficient condition for pattern formation
        if isempty(G)
            conditions(jj,kk) = 0;
        else
            conditions(jj,kk) = max(alpha(kk)*G.*(1-(1-G)/r(kk)) - alpha(jj)*(1-G).*(1-(1-G)/r(kk))...
                - phi(f_0, f_1, G, t_2, s_2) + phi_prime(f_0, f_1, G, t_2, s_2).*(1-G) );
        end
    end
end

%temp = abs(conditions) < 5e-3;
%conditions(temp) = 0;

figure(3);
imagesc(alpha,r,0.5*(sign(conditions')+1));
colormap(bluewhitered);
xlabel('\alpha');
ylabel('r');
set(gca,'YDir','normal');
set(gca,'linewidth',2);
set(gca,'FontSize',36);
toc;