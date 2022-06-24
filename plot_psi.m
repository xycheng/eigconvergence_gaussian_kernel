function plot_psi(fig_num, idx_mu, tt1, psi1, tt2, psi2)


figure(fig_num),clf;
%
subplot(321), hold on;
v1=psi1(:,idx_mu(1));
v2=psi2(:,idx_mu(1));
tmp = interp1(tt2,v2,tt1 ,'spline');
v1 = v1*sign(v1'*tmp);
plot(tt1,v1 ,'-')
plot(tt2,v2 ,'.')
grid on;
title(sprintf('k=%d',idx_mu(1)))

subplot(322), hold on;
plot( tt1, abs(v1-tmp), '.-');
grid on;

%
subplot(323), hold on;
v1=psi1(:,idx_mu(2));
v2=psi2(:,idx_mu(2));
tmp = interp1(tt2,v2,tt1 ,'spline');
v1 = v1*sign(v1'*tmp);
plot(tt1,v1 ,'-')
plot(tt2,v2 ,'.')
grid on;
title(sprintf('k=%d',idx_mu(2)))

subplot(324), hold on;
plot( tt1, abs(v1-tmp), '.-');
grid on;

%
subplot(325), hold on;
v1=psi1(:,idx_mu(3));
v2=psi2(:,idx_mu(3));
tmp = interp1(tt2,v2,tt1 ,'spline');
v1 = v1*sign(v1'*tmp);
plot(tt1,v1 ,'-')
plot(tt2,v2 ,'.')
grid on;
title(sprintf('k=%d',idx_mu(3)))

subplot(326), hold on;
plot( tt1, abs(v1-tmp), '.-');
grid on;

return;