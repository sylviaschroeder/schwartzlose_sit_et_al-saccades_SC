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

saccadeIntervals = [];
amplitudes = [];
vel_stat = [];
onsetXY = [];
if all(isnan(x)) || all(isnan(y))
    return
end

maxStep = 10;

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

% determine saccade onsets and offsets for step between eye positions from 1
% to maxStep samples appart
saccOnsets = cell(maxStep,1);
saccOffsets = cell(maxStep,1);
velocity = NaN(length(x), maxStep);
gaussPars = NaN(2, maxStep);
for step = 1:maxStep
    diffX = padarray(x(1+step : end) - x(1 : end-step), step, 0, 'pre');
    diffY = padarray(y(1+step : end) - y(1 : end-step), step, 0, 'pre');

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

    velocity(:,step) = sqrt(diffX.^2 + diffY.^2);
    % determine suitable velocity threshold to detect saccades: fit Gaussian to
    % distribution of velocity (on  log-log scale), use x STDs (default: 1)
    % as threshold
    vel_log = log(velocity(velocity(:,step)>0, step));

    % NEW CODE TO FIND SACCADE THRESHOLD--------------------------------------
    f = mle(vel_log); % fitting a normal distribution to the distribution of
    % logarithmic velocity values
    gaussPars(:,step) = f;
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

    % detect saccade onsets (with minimum distance of minDist)
    sacc = velocity(:,step) > thresh_vel;
    % set values after each saccade to 1 -> values between nearby saccades 
    % will be >0
    sacc = conv(double(sacc), ones(1, minDist));
    sacc(length(x)+1:end) = [];
    onsets = find(diff(sacc>0) > 0) + 1;
    if step == 1
        onsets = onsets - 1;
        % if velocity one time point before detected threshold crossing is still
        % larger than half the threshold, shift saccade onset back
        velPrev = velocity(onsets-1,step) > thresh_vel/2;
        onsets(velPrev) = onsets(velPrev) - 1;
    end

    % detect saccade offsets
    saccAfter = velocity(:,step) > thresh_vel;
    saccAfter = conv(double(saccAfter), ones(1, minDist));
    saccBefore = velocity(:,step) > thresh_vel;
    saccBefore = conv(double(saccBefore), ones(1, minDist), 'valid');
    saccBefore = [saccBefore; zeros(length(saccAfter) - length(saccBefore),1)];
    saccBetweenClose = saccAfter>0 & saccBefore>0;
    sacc = double(velocity(:,step) > thresh_vel/2 | saccBetweenClose(1:length(velocity)));
    offsets = find(diff(sacc) < 0);
    offsets = offsets' - onsets;
    indOffsets = offsets > 0;
    indOffsets = indOffsets & [true(size(onsets)), ~indOffsets(:,1:end-1)];

    if sum(indOffsets(end, :)) <1 % if the last offset is not detected, remove the last onset
        onsets = onsets(1:end-1);
    end
    offsets = onsets + offsets(indOffsets);

    saccOnsets{step} = onsets;
    saccOffsets{step} = offsets;
end

% determine optimal step size for saccade detection
numSacc = cellfun(@length, saccOnsets);
dNumSacc = [0; diff(numSacc)];
[m, optimum] = max(dNumSacc);
if m < 0.1 * numSacc(1)
    optimum = 1;
end

% if optimal time difference to detect saccades is >1:
% find optimal onset time of each saccade: take 1st time point before 
% detected onset, where velocity between neighbouring time samples is >0 +
% find optimal offset time of each saccade: take last time point before 
% detected offset, where velocity between neighbouring time samples is >0
saccadeIntervals = [saccOnsets{optimum} saccOffsets{optimum}];
if optimum > 1
    lowVel = exp(gaussPars(1,optimum) - gaussPars(2,optimum));
    for s = 1:size(saccadeIntervals,1)
        si = saccadeIntervals(s,:);
        % between time samples used to determine velocity of detected 
        % saccade onset: include first maximum of 1-diff velocity into 
        % saccade
        saccVel1 = velocity(si(1) + (-optimum:0), 1);
        [~,indMax] = max(saccVel1(2:end));
        indMaxGlobal = si(1) - optimum + indMax;
        % find optimal onset
        saccVel1 = saccVel1(1:indMax+1);
        % saccVel1 = velocity(saccadeIntervals(s,1)-optimum : ...
        %     min(saccadeIntervals(s,1), saccadeIntervals(s,2)), 1);
        ind = find(saccVel1 < lowVel,1,'last');
        if isempty(ind)
            ind = 1;
        end
        saccadeIntervals(s,1) = si(1) - (optimum+1) + ind;
        % find optimal offset
        indOffset = max(indMaxGlobal, si(2)-optimum) : si(2);
        saccVel1 = velocity(indOffset, 1);
        ind = find(saccVel1 > lowVel,1,'last');
        if isempty(ind)
            ind = length(indOffset);
        end
        saccadeIntervals(s,2) = indOffset(ind);
    end
end

% check that offset is still more temporal/nasal compared to onset location
diffs = x(saccadeIntervals(:,2)) - x(saccadeIntervals(:,1));
valid = true(size(diffs));
switch dir
    case 'temp'
        valid = diffs < 0;
    case 'nas'
        valid = diffs > 0;
end

saccadeIntervals = saccadeIntervals(valid,:);
onsetXY = [x(saccadeIntervals(:,1)), y(saccadeIntervals(:,1))];
velocity = velocity(:,optimum);
vel_stat.velocity = velocity;
vel_stat.gauss_fit = gaussPars(:,optimum);
thresh_vel = exp(gaussPars(1,optimum) + scaleThresh*gaussPars(2,optimum));

% determine amplitude of each saccade (in pixels)
% amplitudes = zeros(size(onsets));
for am = 1:size(saccadeIntervals,1)
    % calculate pairwaise distance between each sample during saccade
    xSacc = x(saccadeIntervals(am,1) : saccadeIntervals(am,2));
    ySacc = y(saccadeIntervals(am,1) : saccadeIntervals(am,2));
    dists = sqrt((xSacc - xSacc').^2 + (ySacc - ySacc').^2);
    dists_x = xSacc - xSacc';
    dists_y = ySacc - ySacc';

    [amplitudes.vec(am), idx] = max(dists(:));
    amplitudes.x(am) = dists_x(idx);
    amplitudes.y(am) = dists_y(idx);
end

if doPlot > 0
    xl = [1 length(x)];
    figure('WindowState', 'maximized')
    ax = zeros(1,3);
    subplot(3,6,1:5) % x-position
    plot(x,'k')
    hold on
    plot(saccadeIntervals(:,1), x(saccadeIntervals(:,1)), '>r')
    plot(saccadeIntervals(:,2), x(saccadeIntervals(:,2)), '<b')
    ylabel('eye pos in x')
    ax(1) = gca;
    subplot(3,6,7:11) % y-position
    plot(y,'k')
    hold on
    plot(saccadeIntervals(:,1), y(saccadeIntervals(:,1)), '>r')
    plot(saccadeIntervals(:,2), y(saccadeIntervals(:,2)), '<b')
    ylabel('eye pos in y')
    ax(2) = gca;
    subplot(3,6,13:17) % velocity
    plot(vel_stat.velocity,'k')
    hold on
    plot(xl, [1 1].*thresh_vel, 'r')
    plot(saccadeIntervals(:,1), velocity(saccadeIntervals(:,1)), '>r')
    plot(saccadeIntervals(:,2), velocity(saccadeIntervals(:,2)), '<b')
    ylabel('velocity')
    ax(3) = gca;
    linkaxes(ax, 'x')
    xlim(xl)
    xlabel('Time (in samples)')
    
    subplot(3,6,18)
    vel_log = log(velocity(velocity>0));
    histogram(vel_log, 'Normalization', 'pdf')
    hold on
    v = min(vel_log):0.01:max(vel_log);
    plot(v, normpdf(v, gaussPars(1,optimum), gaussPars(2,optimum)), 'r')
    plot([1 1].*log(thresh_vel), ylim, 'r')
    xlim(v([1 end]))
    xlabel('log(velocity)')
    legend off

    sgtitle(sprintf('%s saccades (%d-sample difference)', str, optimum))
end