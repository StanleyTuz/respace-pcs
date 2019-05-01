function [params] = pcs_fit(X)
%PCS_FIT - Fit a periodic cubic spline (pcs) to a set of 1D points.
%Can be used to fit splines to plane or space curves by applying to
% each coordinate separately.
%
% Syntax:  [N] = pcs_fit(X)
%
% Inputs:
%    X - n-by-1 array of 2D coordinates representing curve. Final point
%    should not be a repeat of the first.
%
% Outputs:
%    params - n-by-4 array of pcs parameters, where the jth row are the
%    parameters for the jth subinterval.
%
% Example: 
%    X = cos(0:.1:2*pi)';
%    px = pcs_fit(X);
%
% Dependencies: 
%
% See also: pcs_eval.m, pcs_respace.m, pcs_integrate_al.m
% 
% Stanley Tuznik
% stanley.tuznik@gmail.com
% Apr 2019; Last revision: 30-Apr-2019

if size(X,1) == 1
    X = X';
end

%% Fit the parameters
n = size(X,1);
Ap = diag(ones(1,n-1), 1) + diag(ones(1,n-1), -1) + 4*diag(ones(1,n), 0) + diag(1, n-1) + diag(1, -n+1);
z = 3*(circshift(X, -1) - 2*X + circshift(X, 1));

[c,~] = bicg(Ap, z);

b = (circshift(X, -1) -  X - (2*c + circshift(c, -1))/3)';
d = ((circshift(c, -1) - c)/3)';
c = c';
a = X';

params = [a ; b ; c ; d];