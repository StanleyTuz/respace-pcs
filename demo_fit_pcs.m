
% demo_fit_pcs.m - demonstrate how to fit curves with a periodic cubic
% spline
%

%% Fit a periodic 1D curve
close all
% construct a curve to fit a spline to
t = [0:.1:2*pi]';

x = cos(2*t) .* sin(t);

% fit spline parameters
px = pcs_fit(x);

% new points to evaluate at
t_n = [0:.01:length(x)]';
x_n = pcs_eval(t_n, px);

figure()
plot(t,x,'b.')
hold on
plot(t_n * 2*pi / length(t) + t(1), x_n, 'r-')
plot(t_n * 2*pi / length(t) + t(1) + 2*pi, x_n, 'r-')
legend('original points', 'spline')
title('1D points, spline, and periodicity')
grid()


%% Fit a 2D periodic curve (closed)
close all
%
t = [0:.1:2*pi]';

x = cos(3*t);
y = .2*sin(t);

% fit spline parameters
px = pcs_fit(x);
py = pcs_fit(y);

% new points to evaluate at
t_n = [0:.1:length(x)]';
x_n = pcs_eval(t_n, px);
y_n = pcs_eval(t_n, py);

figure()
plot(x,y,'b.')
axis([min(x)-.2 max(x)+.2 min(y)-.2 max(y)+.2])
hold on
plot(x_n, y_n, 'r-')
legend('original points', 'spline')
title('2D points and spline')
grid()






%% Fit a 3D periodic curve (closed)
close all
%
t = [0:.1:2*pi]';

x = cos(3*t);
y = .2*sin(t);
z = sin(t).*cos(t);

% fit spline parameters
px = pcs_fit(x);
py = pcs_fit(y);
pz = pcs_fit(z);

% new points to evaluate at
t_n = [0:.1:length(x)]';
x_n = pcs_eval(t_n, px);
y_n = pcs_eval(t_n, py);
z_n = pcs_eval(t_n, pz);

figure()
plot3(x,y,z,'b.')
axis([min(x)-.2 max(x)+.2 min(y)-.2 max(y)+.2 min(z)-.2 max(z)+.2])
hold on
plot3(x_n, y_n, z_n, 'r-')
legend('original points', 'spline')
title('2D points and spline')
grid()
