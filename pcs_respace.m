function [X_resp] = pcs_respace(X, N_new)
%PCS_RESPACE - Resample points along an ordered set of planar points such
%that the arc length between the new points is approximately constant.
% Respacing is done by first determining the arc length distances required
% for each new point, then finding corresponding time values for these
% distances by binary searching.
%
% Syntax:  [X_resp] = pcs_respace(X, N_new)
%
% Inputs:
%    X - N-by-2 array of plane coordinates of curve points along which to
%    resample.
%    N_new - number of new points to be placed along curve.
%
% Outputs:
%    X_resp - N_new-by-2 array of respaced plane coordinated.
%
% Example: 
%    X = [cos(0:.1:2*pi); sin(pi/4:.1:2*pi+pi/4)]';
%    X_resp = pcs_respace(X, 100);
%
% Dependencies: pcs_integrate_al.m
%
% See also: pcs_fit.m, pcs_eval.m
% 
% Stanley Tuznik
% stanley.tuznik@gmail.com
% Apr 2019; Last revision: 29-Apr-2019

% Respace the points 
toler = 1e-5;
% Dependencies: pcs_eval, pcs_fit, pcs_integrate

px = pcs_fit(X(:,1));
py = pcs_fit(X(:,2));
N_old = size(X,1);

L = pcs_integrate_al(N_old, px, py); % total arc length

ds = L / N_new; % new point arc length spacing

S = [0:ds:L]; % new arc length points

T = zeros(size(S));

for s_ind = 1:length(S)
    sj = S(s_ind);
    
    % initial enpoints
    ta = 0; % t_1 = 0
    tb = N_old; % t_1' = N
    tc = (ta+tb)/2;
    
    % start loop
    while 1
        tc = (ta+tb)/2;

        fta = pcs_integrate_al(ta, px, py);
        ftb = pcs_integrate_al(tb, px, py);
        ftc = pcs_integrate_al(tc, px, py);

        if abs(ftc - sj) < toler
            break
        elseif abs(fta - sj) < toler
            tc = ta;
            break
        elseif abs(ftb - sj) < toler
            tc = tb;
            break
        end
        
        % search interval update
        if (fta <= sj) && (sj < ftc)
            tb = tc;
        else
            ta = tc;
        end
 
    end
    
    T(s_ind) = tc;
end

vx = pcs_eval(T', px);
vy = pcs_eval(T', py);

% last point is a repeat of the first; discard it
X_resp = [vx(1:end-1,1), vy(1:end-1,1)];




