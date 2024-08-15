function [saccadeIntervals, amplitudes, vel_stat, onsetXY] = ...
    findSaccades(x, y, minDist, scaleThresh, dir, doPlot)
%FINDSACCADES Detects saccade times and amplitudes.
%   [saccadeIntervals, amplitudes,  vel_stat, onsetXY] = FINDSACCADES(x, y, minDist, ...
%    scaleThresh, doPlot) returns on- and offset times of saccades and
%    their amplitudes.
%
%   saccadeIntervals    [saccades x 2]; start and end of each saccade
%                       (given in number of samples from start of x and y 
%                       traces.
%   amplitudes          [saccades x 1]; amplitude of each saccade, defined
%                       as maximum Eucleadian distance between any eye
%                       positions during the saccade.
%  vel_stats
%
%  onsetXY
%
%   x                   [t x 1]; x position of pupil
%   y                   [t x 1]; y position of pupil
%   minDist             double; minimum time interval between two separate
%                       saccades, if two saccades are detected with less
%                       than this interval the second one is treated as
%                       continuation of the first.
%   saccThresh          double; determines the threshold of pupil velocity
%                       used to detect saccades; given in standard
%                       deviations of the log-log distribution of velocity;
%                       optional; default is 1.
%   dir                 string, either {'all, 'temp', 'nas'}
%                       'all' to analyse every direction; 
%                       'temp' to analyse temporal saccades
%                       'nas' to analyse nasal saccades
%                       optional; default is 'all';
%   doPlot              0 or 1; plots detected saccades if 1;
%                       optional; default is 0

if nargin < 4
    scaleThresh = 1;
end

if nargin < 5
    dir = all;
end

if nargin < 6
    doPlot = 0;
end

% get rid of small wiggles and steep "resets" immediately after a saccade
x = medfilt1(x, 3);
y = medfilt1(y, 3);

%determine velocity
diffX = diff(x);
diffY = diff(y);

switch dir
    case 'all'
    case 'temp'
        % set positive x-movements (in nasal direction) to zero
        nas_idx = diffX > 0; 
        diffX(nas_idx) = 0; 
        diffY(nas_idx) = 0;
        str = 'Temporal';
    case 'nas'
        % set negative x-movements (in nasal direction) to zero
        temp_idx = diffX < 0; 
        diffX(temp_idx) = 0; 
        diffY(temp_idx) = 0;
        str = 'Nasal';
end

velocity = sqrt(diffX.^2 + diffY.^2);
% determine suitable velocity threshold to detect saccades: fit Gaussian to
% distribution of velocity (on  log-log scale), use x STDs (default: 1)
% as threshold
vel_log = log(velocity(velocity>0));

% NEW CODE TO FIND SACCADE THRESHOLD--------------------------------------
f = mle(vel_log); % fitting a normal distribution to the distribution of
                  % logarithmic velocity values
thresh_vel = exp(f(1) + scaleThresh*f(2));
% ------------------------------------------------------------------------

% OLD CODE TO FIND SACCADE THRESHOLD--------------------------------------
% edges_vel = floor(min(vel_log)*10)/10 : 0.1 : ceil(max(vel_log)*10)/10; % standardise? 
% bins_vel = edges_vel(1:end-1) + 0.05;
% n = histcounts(vel_log, edges_vel);
% n = log(n);
% ind = ~isinf(n);
% f = fit(bins_vel(ind)', n(ind)', 'gauss1');
% coeffs = coeffvalues(f);
% NOTE: 3rd coefficient of gauss1 model is NOT sigma, it is sqrt(2)*sigma
% thresh_vel = exp(coeffs(2) + coeffs(3)/sqrt(2)) * scaleThresh; % threshold = mean + STD of velocity distribution (in log-log histogram)
% thresh_vel = exp(sum(coeffs(2:3))) * scaleThresh; % threshold = mean + STD of velocity distribution (in log-log histogram)
% ------------------------------------------------------------------------

vel_stat.velocity = velocity;
vel_stat.gauss_fit = f;

% detect saccade onsets (with minimum distance of minDist)
sacc = velocity > thresh_vel;
sacc = conv(double(sacc), ones(1, minDist));
sacc(length(x)+1:end) = [];
% set values after each saccade to 1 -> values betwee nearby saccades will
% be >0
onsets = find(diff(sacc > 0) > 0) + 1;
% if velocity one time point before detected threshold crossing is still
% larger than half the threshold, shift saccade onset back
velPrev = velocity(onsets-1) > thresh_vel/2;
onsets(velPrev) = onsets(velPrev) - 1;

% detect saccade offsets
saccAfter = velocity > thresh_vel;
saccAfter = conv(double(saccAfter), ones(1, minDist));
saccBefore = velocity > thresh_vel;
saccBefore = conv(double(saccBefore), ones(1, minDist), 'valid'); 
saccBefore = [saccBefore; zeros(length(saccAfter) - length(saccBefore),1)];
saccBetweenClose = saccAfter>0 & saccBefore>0;
sacc = double(velocity > thresh_vel/2 | saccBetweenClose(1:length(velocity)));
offsets = find(diff(sacc) < 0) + 1;
offsets = offsets' - onsets;
indOffsets = offsets > 0;
indOffsets = indOffsets & [true(size(onsets)), ~indOffsets(:,1:end-1)];

if sum(indOffsets(end, :)) <1 % if the last offset is not detected, remove the last onset
    onsets = onsets(1:end-1);
end
offsets = onsets + offsets(indOffsets);

saccadeIntervals = [onsets offsets];

onsetXY = [x(onsets), y(onsets)];

% determine amplitude of each saccade (in pixels)
% amplitudes = zeros(size(onsets));
for am = 1:numel(onsets)
    % calculate pairwaise distance between each sample during saccade
    xSacc = x(onsets(am) : offsets(am));
    ySacc = y(onsets(am) : offsets(am));
    dists = sqrt((xSacc - xSacc').^2 + (ySacc - ySacc').^2);
    dists_x = xSacc - xSacc';
    dists_y = ySacc - ySacc';

    [amplitudes.vec(am), idx] = max(dists(:));
    amplitudes.x(am) = dists_x(idx);
    amplitudes.y(am) = dists_y(idx);

end


if doPlot > 0
    t = 1:length(x);
    figure('WindowState', 'maximized')
    ax = zeros(1,3);
    subplot(3,6,1:5) % x-position
    plot(x,'k')
    hold on
    plot(t(onsets), x(onsets), '>r')
    plot(t(offsets), x(offsets), '<b')
    ylabel('eye pos in x')
    ax(1) = gca;
    subplot(3,6,7:11) % y-position
    plot(y,'k')
    hold on
    plot(t(onsets), y(onsets), '>r')
    plot(t(offsets), y(offsets), '<b')
    ylabel('eye pos in y')
    ax(2) = gca;
    subplot(3,6,13:17) % velocity
    plot(2:length(x),velocity,'k')
    hold on
    plot(t([1 end]), [1 1].*thresh_vel, 'r')
    plot(t(onsets), velocity(onsets), '>r')
    plot(t(offsets), velocity(offsets), '<b')
    ylabel('velocity')
    ax(3) = gca;
    linkaxes(ax, 'x')
    xlim(t([1 end]))
    xlabel('Time (in samples)')
    
    subplot(3,6,18)
    histogram(vel_log, 'Normalization', 'pdf')
    hold on
    v = min(vel_log):0.01:max(vel_log);
    plot(v, normpdf(v, f(1), f(2)), 'r')
    plot([1 1].*log(thresh_vel), ylim, 'r')
    xlim(v([1 end]))
    xlabel('log(velocity)')
    legend off

    sgtitle(sprintf('%s saccades', str))
end