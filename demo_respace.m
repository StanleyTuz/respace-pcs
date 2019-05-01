
% demo_respace.m - demonstrate how to compute arc length of a plane curve
% represented by a list of points.
%

%% Compute arc length of a circle of radius 3 (circumference).
close all
% construct the points
t = [0:.1:2*pi]';

x = 3*cos(t);
y = 3*sin(t);

px = pcs_fit(x);
py = pcs_fit(y);

al = pcs_integrate_al(length(x), px, py);
disp(['result: ', num2str(al)])
disp(['exact:  ', num2str(2*pi*3)])
disp(['relative error: ', num2str(norm(al-2*pi*3)/norm(2*pi*3))])

