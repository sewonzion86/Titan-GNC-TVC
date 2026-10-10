function res = Titan_Verify_MC(N,doPlots)
% TITAN_VERIFY_MC  Requirements verification + Monte Carlo for Titan_GNC_TVC.
%   res = Titan_Verify_MC;          % 200 runs, with plots
%   res = Titan_Verify_MC(500,false);
% Needs Titan_GNC_TVC.m on the path. UNTESTED on my side (no MATLAB here) - run N=5 first.

if nargin < 1, N = 200; end
if nargin < 2, doPlots = true; end
d = pi/180;

% ---------------- requirements (limit | verification: every run) ----------
REQ = {
 'R1 Touchdown vertical speed  (-2.0 < vz < 0 m/s)'      @(o,m) o.touchVz > -2.0 && o.touchVz < 0
 'R2 Touchdown lateral speed   (< 1.0 m/s)'               @(o,m) o.touchLatSpeed < 1.0
 'R3 Touchdown tilt            (< 10 deg)'                @(o,m) o.touchTiltDeg < 10
 'R4 Landing error             (< 5 m)'                   @(o,m) o.landErr < 5.0
 'R5 Propellant reserve        (>= 10 kg left)'           @(o,m) o.propLeft >= 10
 'R6 TVC margin                (peak gimbal < 10.8 deg = 90% of 12)' @(o,m) o.maxGimbalDeg < 10.8
 'R7 Terminated by touchdown   (no timeout)'              @(o,m) strcmp(o.status,'TOUCHDOWN')
 'R8 Ignition altitude         (> 100 m)'                 @(o,m) ~isnan(o.ignitionAlt) && o.ignitionAlt > 100
 'R9 Final nav error           (< 5 m position)'          @(o,m) m.navErrEnd < 5.0
};
nR = size(REQ,1);

% ---------------- named stress cases ---------------------------------------
CASES = {
 'nominal',          struct()
 'engine out (engine 2, at 100 m)', struct('engFailIdx',2)
 'dense air (+10%)', struct('rhoScale',1.10)
 'thin air (-10%)',  struct('rhoScale',0.90)
 'wind x3',          struct('windScale',3)
 'Isp -3%',          struct('IspScale',0.97)
 'heavy (+3% mass)', struct('massScale',1.03)
 'TVC misalign 3x',  struct('misalign',3*[0.2;-0.15]*d)
 'gyro bias 3x',     struct('gyroBias',3*[0.010;-0.008;0.006]*d)
 'offset start 150m',struct('x0',[-150;100])
};
fprintf('=== STRESS CASES ===\n');
cs = false(size(CASES,1),1);
for i = 1:size(CASES,1)
    try
        [o,m] = runOne(7,CASES{i,2});
        ok = all(cellfun(@(f) f(o,m),REQ(:,2)));
        cs(i) = ok;
        fprintf('%-34s %s  vz %.2f  lat %.2f  err %.2f m  prop left %.1f kg\n', ...
            CASES{i,1}, tern(ok,'PASS','FAIL'), o.touchVz,o.touchLatSpeed,o.landErr,o.propLeft);
    catch ME
        fprintf('%-34s ERROR: %s\n',CASES{i,1},ME.message);
    end
end

% ---------------- Monte Carlo ------------------------------------------------
fprintf('\n=== MONTE CARLO (N = %d) ===\n',N);
rs = RandStream('mt19937ar','Seed',2026);
M = nan(N,12);  PASS = false(N,nR);  ALL = false(N,1);
for k = 1:N
    u = @(a,b) a + (b-a)*rand(rs);
    g = @() randn(rs);
    op = struct( ...
        'windScale', u(0.5,3.0), ...
        'rhoScale',  u(0.90,1.10), ...
        'massScale', u(0.97,1.03), ...
        'IspScale',  u(0.97,1.03), ...
        'TmaxScale', u(0.90,1.00), ...
        'x0',        [-60;30] + 40*[g();g()], ...
        'v0',        [0.8;0.1;-7] + [0.5*g();0.5*g();0.5*g()], ...
        'tilt0',     abs(3*d + 1.5*d*g()), ...
        'misalign',  [0.2;-0.15]*d + 0.1*d*[g();g()], ...
        'gyroBias',  [0.010;-0.008;0.006]*d .* (1 + 0.5*[g();g();g()]), ...
        'accBias',   [0.004;-0.003;0.005] + 0.002*[g();g();g()]);
    try
        [o,m] = runOne(1000+k,op);
        for r = 1:nR, PASS(k,r) = REQ{r,2}(o,m); end
        ALL(k) = all(PASS(k,:));
        M(k,:) = [o.touchVz o.touchLatSpeed o.touchTiltDeg o.landErr o.propLeft ...
                  o.maxGimbalDeg o.ignitionAlt o.burnTime o.peakQdyn m.navErrEnd o.touchSpeed o.tFlight];
    catch
        % crashed run counts as fail (PASS row stays false)
    end
end

nOK = sum(ALL);  [lo,hi] = wilson(nOK,N);
fprintf('Overall pass: %d/%d = %.1f %%  (95%% Wilson CI %.1f - %.1f %%)\n',nOK,N,100*nOK/N,100*lo,100*hi);
fprintf('\nPer-requirement pass rate:\n');
for r = 1:nR, fprintf('  %-70s %5.1f %%\n',REQ{r,1},100*mean(PASS(:,r))); end

names = {'touchVz (m/s)','lateral (m/s)','tilt (deg)','landErr (m)','propLeft (kg)', ...
         'maxGimbal (deg)','ignitionAlt (m)','burn (s)','peak q (Pa)','navErr (m)','touchSpeed (m/s)','flight (s)'};
fprintf('\n%-18s %9s %9s %9s %9s\n','metric','mean','std','p1','p99');
for j = 1:numel(names)
    x = M(~isnan(M(:,j)),j);
    if isempty(x), continue; end
    fprintf('%-18s %9.3f %9.3f %9.3f %9.3f\n',names{j},mean(x),std(x),prctile1(x,1),prctile1(x,99));
end

res.cases = cs;  res.pass = PASS;  res.all = ALL;  res.metrics = M;  res.metricNames = names;
try
    T = array2table(M,'VariableNames',matlab.lang.makeValidName(names));
    T.pass = ALL;  writetable(T,'titan_mc_results.csv');
    fprintf('\nWrote titan_mc_results.csv\n');
catch
end

if doPlots
    figure('Name','Monte Carlo','Position',[80 80 1300 700]);
    idx = [1 2 3 4 5 6];  lim = [-2 1 10 5 10 10.8];
    for i = 1:6
        subplot(2,3,i); histogram(M(:,idx(i)),30); hold on; grid on;
        xline(lim(i),'r','LineWidth',1.5); title(names{idx(i)});
    end
    figure('Name','Dispersion scatter');
    scatter(M(:,4),-M(:,1),25,double(ALL),'filled'); grid on; colormap([1 0 0; 0 0.6 0]);
    xlabel('Landing error (m)'); ylabel('Descent speed at touchdown (m/s)'); title('Pass (green) / fail (red)');
end
end

% ============================================================================
function [o,m] = runOne(seed,op)
o = Titan_GNC_TVC(false,seed,op);
L = o.log;
m.navErrEnd = norm(L.pE(:,end) - L.X(1:3,end));
end

function s = tern(c,a,b)
if c, s = a; else, s = b; end
end

function [lo,hi] = wilson(k,n)
z = 1.96;  p = k/n;
den = 1 + z^2/n;
c = (p + z^2/(2*n))/den;
h = z*sqrt(p*(1-p)/n + z^2/(4*n^2))/den;
lo = max(0,c-h);  hi = min(1,c+h);
end

function v = prctile1(x,p)
x = sort(x(:));  n = numel(x);
 i = max(1,min(n,round(p/100*n)));
v = x(i);
end
