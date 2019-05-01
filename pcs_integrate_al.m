function [int] = pcs_integrate_al(t, params_x, params_y)
%PCS_INTEGRATE_AL - Compute arc length of a planar pcs up to a specified time.
%Solves the arc length integral for a pcs on [0,t]. Assumes spline knots
%occur at integer values of t.
%The integral is computed over each bin via Gaussian quadrature.
%More curve points (for parameter computation) means closer spacing, which
%in turn means more accurate integration.
%
% Future updates may include arc length calculations for 3D curves.
%
% Syntax:  [N] = pcs_integrate_al(t, params_x, params_y)
%
% Inputs:
%    t - upper limit of arc length integration. Should be on [0,N].
%    params_x, params_y - pcs parameters for spline being integrated.
%
% Outputs:
%    int - scalar arc length value.
%
% Example: 
%    X = [cos(0:.1:2*pi); sin(0:.1:2*pi)]';
%    px = pcs_fit(X(:,1));
%    py = pcs_fit(X(:,2));
%    int = pcs_integrate_al(length(X(:,1)), px, py);
%
% Dependencies: 
%
% See also: pcs_fit.m, pcs_eval.m, pcs_respace.m
% 
% Stanley Tuznik
% stanley.tuznik@gmail.com
% Apr 2019; Last revision: 29-Apr-2019

%t = mod(t, size(params_x,1)); % constrain to a single trace

% which bin is t in?
k = floor(t) + 1;

ki = k - 1; % which bins to integrate over fully
            % If zero, t in [0,1)
            
% how much to integrate in kth bin? Fractional part
t_part = mod(t,1);


% In total, need to integrate over ki bins, then over the fractional part
% t_part in the kth bin.
int = 0;
% integrate over full bins
if (ki ~= 0)
    for kk = 1:ki
        % integrate over the kkth bin, with appropriate params
        % int = int + Integrate(1.0, tkk, px(kk), py(kk))
        int = int + int_bin(1.0, kk-1, params_x(:,kk), params_y(:,kk));
    end
end

% integrate over fractional part
if (t_part ~= 0)
    % int = int + Integrate(t_part, tk, px(k), py(k))
    int = int + int_bin(t_part, k-1, params_x(:,k), params_y(:,k));
end

    
   

function [gc_int] = int_bin(t_frac, tk, pxk, pyk)
% Perform arc length integration via Gaussian quadrature over a single bin.
% Importantly, the function is not piecewise defined over a single bin.
% t_frac is how far along the bin to integrate, in [0,1].
% tk is the time at the start of the bin.
% pxk, pyk are the parameters for the pcs on this bin.
%   pxk = params_x(:,k);
%   pyk = params_y(:,k);
%
% Perform 5-point Gaussian quadrature.

t_frac = tk + t_frac;
f = @(x) sqrt( ( pxk(2) + 2*pxk(3)*(x-tk) + 3*pxk(4)*(x-tk)^2 )^2 +  ...
    ( pyk(2) + 2*pyk(3)*(x-tk) + 3*pyk(4)*(x-tk)^2 )^2 );

% points
u = [ -(1/3)*sqrt(5+2*sqrt(10/7)),
    -(1/3)*sqrt(5-2*sqrt(10/7)),
    0,
    (1/3)*sqrt(5-2*sqrt(10/7)),
    (1/3)*sqrt(5+2*sqrt(10/7))]';

% weights
w = [ (322-13*sqrt(70))/900,
    (322+13*sqrt(70))/900,
    128/225,
    (322+13*sqrt(70))/900,
    (322-13*sqrt(70))/900]';

int = 0;
for j = 1:5 % sum it all up
    int = int + w(j)*f( (t_frac-tk)/2 * u(j) + (t_frac+tk)/2 );
end

gc_int = int * (t_frac-tk)/2;
