tic
%% Numerical method parameters
L = 1; % working on [0,L]
N = 600; % N+1 grid points
delta = L/N;  % spatial discretization parameter
error_tolerance = 1e-14; % for the nonlinear solver
eig_eps = 1e-7; % error tolerance for removing zero eigenvalue
het_tol = 1e-4; % tolerance to classify solution as homo/heterogeneous
display_sol = 1;
%% Function definitions
P_fun = @(x, p_0, p_1) p_0 + p_1.*x;
alpha_p = @(alpha, alpha_s, p) alpha + alpha_s.*p;
J_F_fun = @(x, a, sigma_F) exp( -( (a-x).^2)/(2*sigma_F^2) )/(sqrt(2*pi*sigma_F^2)); %1./(pi*sigma_F*(1+((x-a).^2)/sigma_F^2));%
W_fun = @(x, a, sigma_W) exp( -( (a-x).^2 )/(2*sigma_W^2) )/(sqrt(2*pi*sigma_W^2)); %1./(pi*sigma_W*(1+((x-a).^2)/sigma_W^2));%
phi = @(f_0_fun, f_1, g, t_2_fun, s_2) f_0_fun + (f_1-f_0_fun)./(1 + exp(-(g-t_2_fun)/s_2));
%% Model parameters (values from PNAS paper)
p_left=0;
p_right=0;
p_0=p_left;
p_1=(p_right-p_left)/L;

alpha = 5.6;
alpha_s = 0;

r = 0.84;

f_0_ref = 0.1;
f_1_ref = 0.9;
t_2_ref = 0.4;
s_2 = 0.05; % s_2's standard value is 0.05

% Dispersal parameters
sigma_F = 0.1; % seed dispersal radius forest trees
sigma_W = 0.025; % fire spread radius
sigma_R = 0.15;

disp = 4.8:0.01:4.99;
%disp = 3.6125;%linspace(3.5,3.62,20);
disp = fliplr(disp);
branches = 1;
n_eigs_saved=10;

eig_removed = zeros(branches,length(disp));
sol_norm = zeros(branches,length(disp));
sol_norm_min = zeros(branches,length(disp));
sol_norm_max = zeros(branches,length(disp));
Eigenvalues_Saved= zeros(branches,length(disp),n_eigs_saved);

error = zeros(branches,length(disp));
principal_eig = zeros(branches,length(disp));
fig_count = 0;
%G = load('unstable_branch.mat'); % fill in unstable branch of the saddle
%G = G.G;
for jj = 1:branches % 1 for high grass IC, 2 for low grass IC,
    % 4 for mixed spatial stripes between the two species, 3 for a sigmoidal transition
    for count = 1:length(disp)
        fprintf(['progress = ',num2str((    count-1)/(length(disp))*100),'%']);
        fprintf('\n');
        fig_count = fig_count + 1;
        alpha = disp(count);
        %% Set up the initial distributions of the cover types on the grid
        if count == 1
            %G = 0.445*ones(1,N+1);% + 0.02*(sin(4*pi*(0:delta:L)));%+ 0.02*(cos(5*pi*(0:delta:L)));
            %G = G + 0.01;%*(1+cos(3*pi*(0:delta:L)));
            %G = G(end,N+1:2*N+1); % if using time stepping result as initial condition
            %G = tempGsave;
        end
        %if count == 1
        %0.1.*(1+cos(2*pi*(0:delta:L))) + 0.1.*(1+sin(pi*(0:delta:L)))...
        %+ 0.1.*(1+cos(5*pi*(0:delta:L))) + 0.1.*(1+sin(10*pi*(0:delta:L)))...
        %G = 0.2 + 0.1.*(1+sin(7*2*pi*(0:delta:L))); % sinusoidal initial condition
        %end
        %    +0.05.*(1+cos(5*pi*(0:delta:L))) + 0.05.*(1+sin(10*pi*(0:delta:L))); % sinusoidal initial condition
        %         if jj == 19 && count == 1
        %             G = 1-0.15*(0:delta:L);%(1+cos(12*pi*(0:delta:L)));
        %         elseif jj == 71 && count == 1
        %             G = 0.05*(1+cos(2*pi*(0:delta:L)));
        %         elseif jj == 11 && count == 1
        %             G = 0.15*(1 + cos(4*pi*(0:delta:L)) + sin(6*pi*(0:delta:L)) );
        %         elseif jj == 17 && count == 1
        %             G = phi(0.95, 0.05, 0:delta:L, 0.5, 0.05); % front pinned branch first
        %         elseif jj == 3 && count == 1
        %             G = phi(0.95, 0.6, 0:delta:L, 0.55, 0.15); % trying to find unstable branch connecting saddles
        %         elseif jj == 1 && count == 1
        %             G = 0.44*(0:delta:L).^2 - 0.62*(0:delta:L)+0.29;%phi(0.4, 0.1, 0:delta:L, 0.2, 0.1);%0.2 - 0.1*(0:delta:L); %.*(1+cos(12*pi*(0:delta:L)));
        %         end
        %% plot the initial condition for this simulation
        if display_sol == 1
            figure(fig_count);
            subplot(1,2,1), plot(0:delta:L,G,'-.k','LineWidth',2);
            hold on;
            xlim([0 L]); ylim([0 1]);
            set(gca,'linewidth',1.25); set(gca,'FontSize',18);
            xlabel('\Omega'); ylabel('Density');
            set(gca,'linewidth',1.25); set(gca,'FontSize',12); hold on;
        end
        %keyboard;
        %% set up the spatial grid and pre calculate the convolution matrices
        X = 0:delta:L;
        X_L = X-L;
        X_L = X_L(1,1:end-1);
        X_R = X+L;
        X_R = X_R(1,2:end);
        temp_normalise = ones(1,N+1);
        % Save tempW matrices to avoid computing them again
        tempW = ones(N+1,3*N+1);
        Trap = ones(1,3*N+1);
        Trap(1,3*N+1)=0.5;
        Trap(1,1)=0.5;
        for i = 1:N+1
            tempW(i,:) =  W_fun([X_L X X_R],(i-1)*delta,sigma_W);
            temp_normalise(1,i) = sum(tempW(i,:))*delta;
        end
        C_W = max(temp_normalise);
        % Compute the birth and mortality matrices as a function of rainfall
        P_grad = P_fun(X,p_0,p_1); % compute the rainfall gradient along the x-axis
        alpha_grad = alpha_p(alpha, alpha_s, P_grad);
        tempF = ones(N+1,3*N+1);
        temp_normalise_F = ones(1,N+1);
        tempR = ones(N+1,3*N+1);
        temp_normalise_R = ones(1,N+1);
        %% pre-calculate the convolution matrices
        for k = 1:N+1
            tempF(k,:) = J_F_fun([X_L X X_R],(k-1)*delta,sigma_F);
            temp_normalise_F(1,k) = sum(tempF(k,:))*delta;
            tempR(k,:) = J_F_fun([X_L X X_R],(k-1)*delta,sigma_R);
            temp_normalise_R(1,k) = sum(tempR(k,:))*delta;
        end
        C_F = max(temp_normalise_F);
        C_R = max(temp_normalise_R);
        %% Solve the nonlinear problem f(x) = 0 for the equilibrium solution
        % periodic boundary conditions
        spatial_grid = 0:delta:L;
        f = @(x) -alpha_grad.*((delta*sum( tempF.*(1 - repmat([x(1:N) x x(2:N+1)],N+1,1)).*Trap, 2))'/C_F).*x.*(1 - 1/r + ((delta*sum( tempR.*(repmat([x(1:N) x x(2:N+1)],N+1,1)).*Trap, 2))'/C_R)/r) + ...
            phi(f_0_ref, f_1_ref, ((delta*sum( tempW.*(repmat([x(1:N) x x(2:N+1)],N+1,1)).*Trap, 2))'/C_W), t_2_ref, s_2).*(1-x);
        %f = @(x) -alpha_grad.*((delta*sum( tempF.*(1 - repmat([x(1:N) x x(2:N+1)],N+1,1)).*Trap, 2))'/C_F).*x.*(1 - 1/r + ((delta*sum( tempR.*(repmat([x(1:N) x x(2:N+1)],N+1,1)).*Trap, 2))'/C_R)/r) + ...
        %    phi(f_0_ref, f_1_ref, ((delta*sum( tempW.*(repmat([x(1:N) x x(2:N+1)],N+1,1)).*Trap, 2))'/C_W), t_2_ref, s_2).*(1-x);
        
        % second version of f is the "open" boundary condition
        %f = @(x) - alpha_grad.*((delta*sum( tempF.*(1 - repmat([0*x(1:N) x 0*x(2:N+1)],N+1,1)).*Trap, 2))'/C_F).*x.*(1 - 1/r + ((delta*sum( tempR.*(repmat([0*x(1:N) x 0*x(2:N+1)],N+1,1)).*Trap, 2))'/C_R)/r) + ...
        %    phi(f_0_ref, f_1_ref, ((delta*sum( tempW.*(repmat([0*x(1:N) x 0*x(2:N+1)],N+1,1)).*Trap, 2))'/C_W), t_2_ref, s_2).*(1-x);
        options = optimoptions('fsolve','Display','iter','OptimalityTolerance',...
            error_tolerance,'StepTolerance',error_tolerance);
        [G,fval,exitflag,output,jacobian] = fsolve(f,G,options); % use previous solution value as initial guess for nonlinear solve
        % store the relative error for analysis later
        error(jj,count) = delta*trapz(abs(fval));
        sol_norm(jj,count) = delta*trapz(G); % calculate L^1 norm of the solution
        sol_norm_min(jj,count) = min(G);
        sol_norm_max(jj,count) = max(G);
        %%
        if display_sol == 1
            figure(fig_count);
            %         % fancy plot can be commented out normally
            %         %area(0:delta:L,ones(1,N+1),'LineWidth',1.5,'FaceColor',[0 0.39 0]); % fancy plot
            %         %hold on; % fancy plot
            %         %area(0:delta:L,G,'LineWidth',1.5,'FaceColor',[0.565 0.933 0.565]); % fancy plot
            subplot(1,2,1), plot(0:delta:L,G,'-.b','LineWidth',2);
            xlim([0 L]); ylim([0 1]);
            xlabel('\Omega'); ylabel('Density');
            %         %legend('Forest','Grass');
            set(gca,'linewidth',1.25); set(gca,'FontSize',18);
            legend('IC','final soln');
            title(['\alpha = ',num2str(disp(count)),' and L^1 error = ',num2str(error(jj,count))]);
        end
        z = eig(jacobian);
        % remove zero eigenvalue due to translation invariance
        temp_eig = real(z( abs(real(z))>eig_eps ));
        eig_removed(jj,count) = sum(size(z) - size(temp_eig)); % should always be 1
        z = temp_eig; % replace with the adjusted eigenvalues
        principal_eig(jj,count) = max(temp_eig);
        %principal_eig(jj,count) = max(real(z));
        [a,b]=sort(abs(real(z)));
        
        Eigenvalues_Saved(jj,count,:)=z(b(1:n_eigs_saved));
        %Eigenvalues_Saved(jj,count,:) = temp_eig(b(1:n_eigs_saved));
        if display_sol == 1
            subplot(1,2,2), scatter(real(z),imag(z));
            xlim([min(real(z))-0.05 max(real(z))+0.05]);
            ylim([min(imag(z))-0.025 max(imag(z))+0.025]);
            title(['Principal eigenvalue: ',num2str(principal_eig(jj,count))]);
            hold on;
            subplot(1,2,2), xline(0,'-.r','LineWidth',2);
        end
    end
end
%% Plot the L^1 norm of the solution versus the dispersal parameter sigma
figure(100);
hold on;
L1_tolerance = 1e-5;
for i = 1:branches
    for j = 1:length(disp)
        if principal_eig(i,j) < 0 && error(i,j) < L1_tolerance % stable, meets solver tolerances
            if abs(sol_norm_max(i,j)-sol_norm_min(i,j)) > het_tol
                scatter(disp(j),sol_norm(i,j),'g','filled');%,"MarkerEdgeColor","k");
                %scatter(disp(j),sol_norm_max(i,j),'g','filled');
                %scatter(disp(j),sol_norm_min(i,j),'g','filled');
            else
                scatter(disp(j),sol_norm_max(i,j),'r','filled');
            end
        elseif principal_eig(i,j) < 0 && error(i,j) > L1_tolerance % stable(?), possible error
            if abs(sol_norm_max(i,j)-sol_norm_min(i,j)) > het_tol
                scatter(disp(j),sol_norm(i,j),'c');
                %scatter(disp(j),sol_norm_max(i,j),'c');
                %scatter(disp(j),sol_norm_min(i,j),'c');
            else
                scatter(disp(j),sol_norm_max(i,j),'m');
            end
        elseif principal_eig(i,j) > 0 && error(i,j) > L1_tolerance % unstable(?), possible error
            if abs(sol_norm_max(i,j)-sol_norm_min(i,j)) > het_tol
                scatter(disp(j),sol_norm(i,j),'c');
                %scatter(disp(j),sol_norm_max(i,j),'c');
                %scatter(disp(j),sol_norm_min(i,j),'c');
            else
                scatter(disp(j),sol_norm_max(i,j),'m');
            end
        else % unstable, meets solver tolerances
            if abs(sol_norm_max(i,j)-sol_norm_min(i,j)) > het_tol
                scatter(disp(j),sol_norm(i,j),'b','filled');
                %scatter(disp(j),sol_norm_max(i,j),'b','filled');
                %scatter(disp(j),sol_norm_min(i,j),'b','filled');
            else
                scatter(disp(j),sol_norm_max(i,j),'k','filled');
            end
        end
    end
end
hold on;
%xlabel('\alpha');
%ylabel('||G||');
xlim([4 9]);
ylim([0.2 0.35]);
set(gca,'linewidth',2);
set(gca,'FontSize',24);
grid on;
toc

%% Leading eigenvalue output

% [a,b]=sort(real(z));
% a(end-10:end);
% Style = {'*','o','.'};
% for eee=1:n_eigs_saved
%     figure(101)
%     hold on
%     plot(disp,real(squeeze(Eigenvalues_Saved(1,:,eee))),Style{mod(eee,length(Style))+1},'MarkerSize',10)
%
% end
% plot(disp,0*real(squeeze(Eigenvalues_Saved(1,:,1))),'k--');

%% Code for saving figures
%figure('codim1_diagram');
%saveas(gcf,'codim1_alpha_highSigmaR_r=0.84_aug.fig');
%saveas(gcf,'codim1_alpha_highSigmaR_r=0.9_jan12.fig');
%openfig('codim1_alpha_highSigmaR_r=0.9.fig')
% yticks([0 0.2 0.4 0.6 0.8 1])
%print('codim1_alpha_oct7th_subcritical','-dpdf','-fillpage') % or use the
% -bestfit option


% dim = [.2 .5 .3 .3];
% str = 'L = 1, r = 0.8, \sigma_{F} = 0.1, \sigma_{W} = 0.025, \sigma_{R} = 0.05';
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
% set(gca,'FontSize',24);

%% Eigenvector code
% [vvv,ddd]=eig(jacobian);
% [a,b]=sort(abs(real(z)));
% %close all
% for i = 3:4
%     figure(102)
%     hold on
%     plot(real(vvv(:,b(i))),'.')
%     figure(103)
%     hold on
%     plot(imag(vvv(:,b(i))),'.')
%     pause()
% end
%
% figure;plot(diff(sum(abs(vvv),1)),'.')
% figure;plot(diff(a),'.')
% figure;plot(real(vvv(:,b(3))),'.')
% figure;plot(real(vvv(:,b(4))),'r.')