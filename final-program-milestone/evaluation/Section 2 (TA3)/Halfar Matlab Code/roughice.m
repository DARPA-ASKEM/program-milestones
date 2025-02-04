% ROUGHICE   Demonstrate that the SIA code siaflat.m can deal adaptively
%            with terrible surface topography.

% space-time grid parameters
J = 60;  K = J;
Lx = 500e3;  Ly = Lx;
secpera = 31556926;
x = linspace(-Lx,Lx,J+1);  y = x;
[xx,yy] = meshgrid(x,y);

% construct and plot the worst possible ice sheet
% you can replace this section with your real data
H0 = 3000 + 1000 * rand(J+1,K+1);
H0(xx.^2 + yy.^2 > 300000^2) = 0;
figure(1),  surf(x/1000,y/1000,H0,'edgealpha',0.4)
xlabel('x  (km)'), ylabel('y  (km)'), zlabel('h  (m)')
title('initial surface')

% run SIA code
tfyears = 50;  dtyears = 0.2;
[H,dtlist] = siaflat(Lx,Ly,J,K,H0,dtyears*secpera,tfyears*secpera);

% show final state
figure(2),  surf(x/1000,y/1000,H,'edgealpha',0.4)
xlabel('x  (km)'), ylabel('y  (km)'), zlabel('h  (m)')
title(['final surface at t_f = ' num2str(tfyears) ' (a)'])

% show adaptive time-stepping
figure(3),  plot(dtlist/secpera,'o')
xlabel('time steps (count)'), ylabel('\Delta t  (a)')