function verifysia(J)
% compare the Halfar (1983) similarity solution
% at year t = 10000 to the numerical solution using
% siaflat(), for a run from t = 200 to t = 20000
secpera = 31556926;
if nargin<1, J=40; end;

L = 1200e3;
dx = 2 * L / J;
[x,y] = meshgrid(-L:dx:L, -L:dx:L);

t1 = 200;
t2 = 20000;
radius=10000;
dist_from_center = sqrt(x.^2 + y.^2);
H0 = sqrt(radius^2 - dist_from_center.^2);
H0(dist_from_center > radius) = 0;
[H2approx,dtlist] = siaflat(L,L,J,J,H0,10.0*secpera,(t2-t1)*secpera);

% side-by-side comparison of numerical and exact result:
figure(2), surf(x/1000,y/1000,H2approx), shading('flat')
xlabel('x (km)'), ylabel('y (km)'), zlabel('numerical thickness (m)')

% figure showing adaptive time-stepping:
figure(3), plot(dtlist / secpera,'o')
xlabel('step'), ylabel('length of step in years')
