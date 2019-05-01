function [vals] = pcs_eval(t_new, params)
%PCS_EVAL - Evaluate a periodic cubic spline at new time values.
%
% Syntax:  [N] = pcs_eval(t_new, params)
%
% Inputs:
%    params - an N-by-4 array of periodic cubic spline parameters where the
%    jth row are the coefficients for the jth subinterval.
%    t_new - a list of time points to evaluate the spline at, assuming the
%    original spline spans the interval [0,N].
%
% Outputs:
%    vals - n-by-1 array of new coordinates.
%
% Example: 
%    px = pcs_fit(x);
%    t = [0:.1:length(x)]';
%    vx = pcs_eval(t, px);
%
% Dependencies: 
%
% See also: pcs_fit.m, pcs_respace.m, pcs_integrate_al.m
% 
% Stanley Tuznik
% stanley.tuznik@gmail.com
% Apr 2019; Last revision: 29-Apr-2019

a = params(1, :)';
b = params(2, :)';
c = params(3, :)';
d = params(4, :)';

np = size(a, 1);

k = mod(floor(t_new), np) + 1;
tk = mod(t_new - k + 1, np);

vals = a(k, 1) + b(k, 1).*tk + c(k, 1).*(tk.^2) + d(k, 1).*(tk.^3);