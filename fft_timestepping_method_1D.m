tic
%% Model Parameters
alpha = 8;
R = 0.6;
f_0 = 0.1;
f_1 = 0.9;
t_2 = 0.4;
s_2 = 0.05;
%% kernel parameters
sigmaJ = 0.1; % seed spread radius
sigmaR = 0.1; % resource competition radius
sigmaW = 0.025;
%% Define the kernels & the sigmoid for fire spread
phi = @(f_0_fun, f_1, g, t_2_fun, s_2) f_0_fun + (f_1-f_0_fun)./(1 + exp(-(g-t_2_fun)/s_2));
kernel = @(x, sigma) (1/(sigma*sqrt(2*pi))).*exp( -(x.^2)/(2*sigma*sigma) );
%% Fourier Points & initial conditions
N = 2^9;
delta = 1; % spatial discretization parameter
X = linspace(-pi,pi,N);
% compute ft's of the kernels on the mesh
kernel_J = fft(kernel(X,sigmaJ)/sum(kernel(X,sigmaJ)));
kernel_R = fft(kernel(X,sigmaR)/sum(kernel(X,sigmaR)));
kernel_W = fft(kernel(X,sigmaW)/sum(kernel(X,sigmaW)));

% initial conditions
G_old = 0.7 + 0.025*rand(1,N);
G_new = zeros(1,N);

%% Explicit Euler in space
h = 0.025;
T = 50000;
figure(2);
for i = 0:T
    progressbar(i,T);
    dU = (1-G_old).*phi(f_0, f_1, delta*ifft(fft(G_old).*kernel_W), t_2, s_2)...
        - alpha*G_old.*(1 - delta*ifft(fft(G_old).*kernel_J)).*(1- 1/R +...
        delta*ifft(fft(G_old).*kernel_R)/R);
    G_new = G_old + h*dU;
    %G_new = min(max(G_new,0),1);
%     if i == 0
%         plot(X,G_old,'k','LineWidth', 2)
%         ylim([0 1]);
%         xlim([min(X) max(X)])
%         hold on
%         drawnow
%     end
    if mod(i,T) == 0 && i > 0 
        plot(X,G_new,'b', 'LineWidth', 3)
        ylim([0 1]);
        xlim([min(X) max(X)])
        hold on
        drawnow
        set(gca,'FontSize',24);
        grid on;
    end
    if mod(i,T) == 0
        fprintf('\n');
        fprintf('Relative error in l^1 norm:\n');
        fprintf(['G: ',num2str(max(abs(G_new-G_old))/max(G_old)),'\n']);
    end
    G_old = G_new;
end
toc


