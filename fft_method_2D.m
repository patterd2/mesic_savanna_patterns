tic

%% Model Parameters
alpha = 6;
R = 1;
f_0 = 0.1;
f_1 = 0.9;
t_2 = 0.4;
s_2 = 0.05;

%% kernel parameters
sigmaJ = 0.1; % seed spread radius
sigmaR = 0.5; % resource competition radius
sigmaW = 0.01;

%% Define the sigmoid for fire spread
phi = @(f_0_fun, f_1, g, t_2_fun, s_2) f_0_fun + (f_1-f_0_fun)./(1 + exp(-(g-t_2_fun)/s_2));
%% Define the kernel
kernel = @(x, y, sigma) (1/(sigma*sigma*2*pi)).*exp( -(x.^2 + y.^2)/(2*sigma*sigma) );
%% Fourier Points
N = 2^9;
X = linspace(-pi,pi,N);
Y = linspace(-pi,pi,N);

% compute ft's of the kernels on the mesh
[x1, y1] = meshgrid(X, Y);
kernel_J = fft2(kernel(x1,y1,sigmaJ)/sum(sum(kernel(x1,y1,sigmaJ))));
kernel_R = fft2(kernel(x1,y1,sigmaR)/sum(sum(kernel(x1,y1,sigmaR))));
kernel_W = fft2(kernel(x1,y1,sigmaW)/sum(sum(kernel(x1,y1,sigmaW))));

% initial conditions
G_old = 0.15*ones(N,N)+0.05*rand(N,N);
G_new = zeros(N,N);

%% Explicit Euler in space
h = 0.05;
T = 10000;
figure(7);
for i = 0:T
    progressbar(i,T);
    dU = (1-G_old).*phi(f_0, f_1, real(ifft2(fft2(G_old).*kernel_W)), t_2, s_2)...
        - alpha*G_old.*(1 - real(ifft2(fft2(G_old).*kernel_J))).*(1-1/R...
        + real(ifft2(fft2(G_old).*real(kernel_R)))/R);
    G_new = G_old + h*dU;
    if mod(i,250) == 0 && i>0
        progressbar(i,T);
        pcolor(G_new)
        colorbar;
        caxis([0 1]);
        shading interp
        colorbar;
        hold on
        drawnow
    end
    if mod(i,T) == 0 && i > 0
        fprintf('\n');
        fprintf('Relative error in l^1 norm:\n');
        fprintf(['G: ',num2str(max(abs(G_new-G_old))/max(G_old)),'\n']);
    end
    G_old = G_new;
end
toc


