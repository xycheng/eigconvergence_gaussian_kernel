% this file computes the error of eigenvalue and eigenctors of graph
% Laplacian matrix built using a gaussian kernel on data ponts on S^1, over
% a range of values of N (number of samples) and epsW (kernel bandwidth).
% It reproduces Fig 2 (and Fig A.1) in the paper
% 
% X. Cheng and N. Wu. "Eigen-convergence of Gaussian kernelized graph 
% Laplacian by manifold heat interpolation". Accepted to ACHA 2022.
% Preprint: arXiv: 2101.09875.

clear all; rng(2021);

%%
if_plot = 0;
    % set if_plot = 1 to generate plots of eigenvalues and eigenvectors of
    % the graph Laplacian, set if_plot = 0 to skip plotting in the loop


%% data and p

p0 = 1; %unit length S1

p0_func = @(t, w, ome) 1+ w*sin( ome * 2*pi *t);
F0_func = @(t, w, ome) t + w*(1- cos(ome * 2*pi * t))/(ome*2*pi);
dp0_func = @(t,w,ome) w*cos( ome * 2*pi *t)*(ome * 2*pi);
d2p0_func = @(t,w,ome) -w*sin( ome * 2*pi *t)*((ome * 2*pi)^2);

% uniform
p_func = @(t) p0_func(t, 0, 1);
F_func = @(t) F0_func(t, 0, 1);
dp_func = @(t) dp0_func(t, 0, 1);
d2p_func = @(t) d2p0_func(t, 0, 1);


%%
dM =1;
omegaM = 3;
map_to_RD_func = @(t) 1/(sqrt(5)*2*pi)*[...
                       cos(2*pi * t), ...
                       sin(2*pi * t), ...
                       2/omegaM*cos(2*pi * omegaM*t), ...
                       2/omegaM*sin(2*pi * omegaM*t)];

%% vis

n_vis = 1000;

t_vis = sort(rand( n_vis,1),'ascend');
data_vis = map_to_RD_func(t_vis);

figure(5),clf;
subplot(121),
plot( t_vis, p_func( t_vis), '.-');
grid on;
subplot(122),
scatter3(data_vis(:,1), data_vis(:,2), data_vis(:,3), 80, p_func(t_vis), 'o', 'filled');
grid on;
colorbar(); view(3);


%% population eigenfunctions and eigenvalues

% q is uniform
q_func = @(t) p0_func(t, 0, 2);
dq_func = @(t) dp0_func(t, 0, 2);
d2q_func = @(t) d2p0_func(t, 0, 2);


maxl = 10;
maxk = maxl*2+1;
ll = zeros(1, maxk);
for l = 1: maxl
    ll(2*l:2*l+1) = l;
end


%% use Gaussian kernel function

tildemh = 1; %m2/(2*m0)
m0h = 1;
m2h = 2;
h_func = @(r) exp(-r/4)/( (4*pi)^(dM/2));


%% grid of N and eps


% grid of N
nrow = 10;
Nx_list = floor(10.^(2.75:0.05:2.75+0.05*(nrow-1)))';


% grid of epsW
ncol = 80;
lb = -4.0; rb = -2.8;

hb = (rb-lb)/ncol;
epsW_list = 10.^(lb+hb:hb:rb)';

nrun = 5;
    % set nrun=5, takes about 700 seconds to finish. Due to less number of
    % replicas, the results have more noise. 
    % set nrun=500, reproduce Fig 2 and Fig A.1 in the paper, takes longer 
    % time to run the code

%%

tic,

err_vals = zeros( maxk, nrow, ncol, nrun);
err_vecs = zeros( maxk, nrow, ncol, nrun);

for irow = 1: nrow
    
    Nx = Nx_list(irow);
    
    fprintf('Nx=%d:', Nx);
    
    for icol=1: ncol
        epsW = epsW_list(icol);
         
        fprintf('epsW=%6.4e, ', epsW);
        
        for irun=1:nrun
            
            %% sample X
            
            tX = sort(rand(Nx,1),'ascend');
            dataX = map_to_RD_func(tX);
          
            %% popoulation eigenfunctions on the data
            psi_nt = zeros(Nx, maxk);
            mu1nt = zeros(1,maxk);
            psi_nt(:,1) =1;
            for l = 1:maxl
                psi_nt(:, ll==l) = [cos(tX*2*pi*l), sin(tX*2*pi*l)]*sqrt(2);
                mu1nt(:, ll==l) = (2*pi*l)^2;
            end
            
            %% graph laplacian
            disXX2 = squareform( pdist(dataX).^2);
            
            W = (epsW^(- dM/2))*h_func( disXX2/epsW);
            dW = sum(W,2);
               
            [v,d]=eig( diag(dW)-W, diag(dW));
            [lam1, tmp]=sort(diag(d),'ascend');
            v1 = v(:,tmp);
            
            %%   normalize
            % lamk is eigenvalue of (I-D^{-1}W), L_rw = 1/eps*(I-D^{-1}W)
            lam1 = lam1(1:maxk);
            lam1= lam1/(epsW*tildemh);
            

            v1 = v1(:,1:maxk);
            v1 = v1*sqrt(p0*Nx);
            
            % phik = rho_X(  psi_k/sqrt(pN) )
            phi1 = psi_nt(:,1:maxk);
            phi1 = phi1 /sqrt(p0*Nx);
            
            
            mu1 = mu1nt(1: maxk)';
            
            %% rotate eigenvectors to align, due to multiplicity 2
            v1(:,1) =v1(:,1)*sign( sum(v1(:,1) ));
            
            for l= 1 : maxl
                idxk = find(ll==l);
                
                u = v1(:,idxk); %Nx-by-2
                uphi = phi1(:,idxk);
                
                rot = u\uphi;
                [uu,ss,vv]=svd(rot);
                rot = uu*vv';
                
                u = u*rot;
                v1(:,idxk) =u;
            end
             %% allow a scalor alphak multiply to v1
            norm_beta = zeros(maxk,1);
            for j=1:maxk
                u = v1(:,j);
                uphi = phi1(:,j);
                norm_beta(j)= norm(uphi,2)/norm(u,2);
                v1(:,j) = v1(:,j)*norm_beta(j);
            end
            
            %% compute relative error
            
            err_val = (lam1(1:maxk)-mu1(1:maxk))./mu1(1:maxk);
            err_val (1) = 0;
            err_vec = (sqrt(sum((v1(:,1:maxk)-phi1(:,1:maxk)).^2,1))./sqrt(sum(phi1(:,1:maxk).^2,1)))';
            
            
            err_vals(:, irow, icol, irun) = err_val;
            err_vecs(:, irow, icol, irun) = err_vec;
            
            
            %% plot
            
            if if_plot && irun==1
                
                idx_mu = [2:maxk ];
                
                figure(1),clf;
                subplot(131), hold on;
                plot(idx_mu, mu1(idx_mu), 'x-');
                plot(idx_mu, lam1(idx_mu), 'o-');
                grid on; xlabel('k')
                subplot(132), hold on;
                plot(idx_mu, abs( mu1(idx_mu)-lam1(idx_mu) )./mu1(idx_mu), 'x-');
                grid on; title('rel error')
                subplot(133), hold on;
                plot( idx_mu(2:end), diff( mu1(idx_mu)), 'x-');
                plot( idx_mu(2:end), diff( lam1(idx_mu)), 'x-');
                grid on; title('eigen gap')
                
                fignum = 2;
                idx_psi_to_plot = [5,6,7];%[3,4,5];
                plot_psi(fignum, idx_psi_to_plot, tX, v1, tX, phi1);
                
                
                
                figure(3),clf; hold on;
                plot(err_val,'x-')
                plot(err_vec,'x-')
                grid on;
                title(sprintf('Nx=%d,sqrt(eps)=%f', Nx, sqrt(epsW)))
                
                drawnow();
            end
            
        end
        
    end
    fprintf('\n');    
end

toc

%% plotting

% eigenvalue
merr = squeeze( mean( mean(abs(err_vals(2:9,:,:,:)), 4),  1)); 

logmerr1 = log10(merr);

% select again
[logminerr1, idx_eps] = min( logmerr1, [],2);
eps_select = epsW_list(idx_eps);
        
figure(20),clf;
idxW = find( log10(epsW_list) < -3.3);

subplot(121), hold on;
imagesc( log10(Nx_list), log10(epsW_list(idxW)), logmerr1(:,idxW)');
plot( log10(Nx_list), log10( eps_select),'x-r');
p = polyfit( log10( Nx_list), log10( eps_select), 1);
plot( log10(Nx_list), log10(Nx_list)*p(1)+p(2),'--g');
 
dim = [0.25 .58 .3 .3];
str = sprintf('slope = %6.4f', p(1));
annotation('textbox',dim,'String',str,'FontSize',20, ...
    'FitBoxToText','on','Interpreter','latex');

xlabel('$\log_{10} N$', 'Interpreter','latex')
ylabel('$\log_{10} \epsilon$', 'Interpreter','latex')
axis([log10(Nx_list(1)), log10(Nx_list(end)), min(log10(epsW_list(idxW))), max(log10(epsW_list(idxW)))])
colorbar();
set(gca,'FontSize',20);
title('$S^1$, $L_{rw}$', 'Interpreter','latex')

subplot(122), hold on;
plot( log10(Nx_list), logminerr1,'x-r'); 
grid on;
p = polyfit( log10( Nx_list), logminerr1, 1);
plot( log10(Nx_list), log10(Nx_list)*p(1)+p(2),'--b');
dim = [0.72 .58 .3 .3];
str = sprintf('slope = %6.4f', p(1));
annotation('textbox',dim,'String',str,'FontSize',20, ...
    'FitBoxToText','on','Interpreter','latex');
xlabel('$\log_{10} N$', 'Interpreter','latex')
set(gca,'XLim', [log10(Nx_list(1)), log10(Nx_list(end))])
set(gca,'FontSize',20);

title('Error of eigenvalue', 'Interpreter','latex')


% eigevector
merr = squeeze( mean( mean(abs(err_vecs(2:9,:,:,:)), 4),  1));
logmerr1 = log10(merr);

% select again
[logminerr1, idx_eps] = min( logmerr1, [],2);
eps_select = epsW_list(idx_eps);
        
figure(21),clf;
idxW = find( log10(epsW_list) > -3.3);

subplot(121), hold on;
imagesc( log10(Nx_list), log10(epsW_list(idxW)), logmerr1(:,idxW)');
plot( log10(Nx_list), log10( eps_select),'x-r');
p = polyfit( log10( Nx_list), log10( eps_select), 1);
plot( log10(Nx_list), log10(Nx_list)*p(1)+p(2),'--g');
 
dim = [0.25 .58 .3 .3];
str = sprintf('slope = %6.4f', p(1));
annotation('textbox',dim,'String',str,'FontSize',20, ...
    'FitBoxToText','on','Interpreter','latex');

xlabel('$\log_{10} N$', 'Interpreter','latex')
ylabel('$\log_{10} \epsilon$', 'Interpreter','latex')
axis([log10(Nx_list(1)), log10(Nx_list(end)), min(log10(epsW_list(idxW))), max(log10(epsW_list(idxW)))])
colorbar();
set(gca,'FontSize',20);

subplot(122), hold on;
plot( log10(Nx_list), logminerr1,'x-r'); 
grid on;
p = polyfit( log10( Nx_list), logminerr1, 1);
plot( log10(Nx_list), log10(Nx_list)*p(1)+p(2),'--b');
dim = [0.72 .58 .3 .3];
str = sprintf('slope = %6.4f', p(1));
annotation('textbox',dim,'String',str,'FontSize',20, ...
    'FitBoxToText','on','Interpreter','latex');
xlabel('$\log_{10} N$', 'Interpreter','latex')
set(gca,'XLim', [log10(Nx_list(1)), log10(Nx_list(end))])
set(gca,'FontSize',20);

title('Error of eigenvector', 'Interpreter','latex')


%% with smooth

merr = squeeze( mean( mean(abs(err_vals(2:9,:,:,:)), 4),  1)); 

logmerr1 = zeros( size(merr));

smooth_trunc = 12; %16;
kk1 = ifftshift(-ncol:ncol-1);

for ii=1:nrow
    x = log10(merr(ii,:));
    x = [x, x(end:-1:1)];
    hx = fft(x);
    hx( abs(kk1) > smooth_trunc) = 0;
    x1 =ifft( hx,'symmetric');
    logmerr1(ii,:)=x1(1:ncol);
end

% select again
[logminerr1, idx_eps] = min( logmerr1, [],2);
eps_select = epsW_list(idx_eps);
        
figure(22),clf;
idxW = find( log10(epsW_list) < -3.3);

subplot(121), hold on;
imagesc( log10(Nx_list), log10(epsW_list(idxW)), logmerr1(:,idxW)');
plot( log10(Nx_list), log10( eps_select),'x-r');
p = polyfit( log10( Nx_list), log10( eps_select), 1);
plot( log10(Nx_list), log10(Nx_list)*p(1)+p(2),'--g');
 
dim = [0.25 .58 .3 .3];
str = sprintf('slope = %6.4f', p(1));
annotation('textbox',dim,'String',str,'FontSize',20, ...
    'FitBoxToText','on','Interpreter','latex');

xlabel('$\log_{10} N$', 'Interpreter','latex')
ylabel('$\log_{10} \epsilon$', 'Interpreter','latex')
axis([log10(Nx_list(1)), log10(Nx_list(end)), min(log10(epsW_list(idxW))), max(log10(epsW_list(idxW)))])
colorbar();
set(gca,'FontSize',20);
title('$S^1$, $L_{rw}$', 'Interpreter','latex')

subplot(122), hold on;
plot( log10(Nx_list), logminerr1,'x-r'); 
grid on;
p = polyfit( log10( Nx_list), logminerr1, 1);
plot( log10(Nx_list), log10(Nx_list)*p(1)+p(2),'--b');
dim = [0.72 .58 .3 .3];
str = sprintf('slope = %6.4f', p(1));
annotation('textbox',dim,'String',str,'FontSize',20, ...
    'FitBoxToText','on','Interpreter','latex');
xlabel('$\log_{10} N$', 'Interpreter','latex')
set(gca,'XLim', [log10(Nx_list(1)), log10(Nx_list(end))])
set(gca,'FontSize',20);

title('Error of eigenvalue', 'Interpreter','latex')


% eigevector
merr = squeeze( mean( mean(abs(err_vecs(2:9,:,:,:)), 4),  1));
logmerr1 = zeros( size(merr));

kk1 = ifftshift(-ncol:ncol-1);

for ii=1:nrow
    x = log10(merr(ii,:));
    x = [x, x(end:-1:1)];
    hx = fft(x);
    hx( abs(kk1) > smooth_trunc) = 0;
    x1 =ifft( hx,'symmetric');
    logmerr1(ii,:)=x1(1:ncol);
end

% select again
[logminerr1, idx_eps] = min( logmerr1, [],2);
eps_select = epsW_list(idx_eps);
        
figure(23),clf;
idxW = find( log10(epsW_list) > -3.3);

subplot(121), hold on;
imagesc( log10(Nx_list), log10(epsW_list(idxW)), logmerr1(:,idxW)');
plot( log10(Nx_list), log10( eps_select),'x-r');
p = polyfit( log10( Nx_list), log10( eps_select), 1);
plot( log10(Nx_list), log10(Nx_list)*p(1)+p(2),'--g');
 
dim = [0.25 .58 .3 .3];
str = sprintf('slope = %6.4f', p(1));
annotation('textbox',dim,'String',str,'FontSize',20, ...
    'FitBoxToText','on','Interpreter','latex');
xlabel('$\log_{10} N$', 'Interpreter','latex')
ylabel('$\log_{10} \epsilon$', 'Interpreter','latex')
axis([log10(Nx_list(1)), log10(Nx_list(end)), min(log10(epsW_list(idxW))), max(log10(epsW_list(idxW)))])
colorbar();
set(gca,'FontSize',20);

subplot(122), hold on;
plot( log10(Nx_list), logminerr1,'x-r'); 
grid on;
p = polyfit( log10( Nx_list), logminerr1, 1);
plot( log10(Nx_list), log10(Nx_list)*p(1)+p(2),'--b');
dim = [0.72 .58 .3 .3];
str = sprintf('slope = %6.4f', p(1));
annotation('textbox',dim,'String',str,'FontSize',20, ...
    'FitBoxToText','on','Interpreter','latex');
xlabel('$\log_{10} N$', 'Interpreter','latex')
set(gca,'XLim', [log10(Nx_list(1)), log10(Nx_list(end))])
set(gca,'FontSize',20);

title('Error of eigenvector', 'Interpreter','latex')


%%
return;
