function D = repDetector_x4_strongChannel(X, Fs, opts)
% repDetector_x4_strongChannel
%
% Reference repetition detector for CW2 (5CCE2SAS – Signals and Systems)
%
% Detects phonation bursts using an energy-based method on the stronger
% microphone channel and returns burst start/end indices and fixed-length
% analysis windows.
%
% Inputs
%   X   : [N x 2] signal (two microphone channels)
%   Fs  : sampling frequency (e.g., 20000)
%   opts (optional struct with fields):
%       thrFrac   (default 0.02)  - threshold fraction of peak activity (x.^4)
%       smoothMs  (default 10)    - smoothing window for activity (ms)
%       minOnMs   (default 30)    - minimum burst duration to keep (ms)
%       minOffMs  (default 20)    - minimum silence gap; shorter gaps are filled (ms)
%       winSec    (default 0.5)   - analysis window length (seconds)
%       maxReps   (default 4)     - keep up to this many bursts (longest)
%
% Output struct D
%   D.chUsed        : 1 or 2 (channel used for detection)
%   D.starts        : burst start indices (samples)
%   D.ends          : burst end indices (samples)
%   D.midpoints     : midpoint of each burst (samples)
%   D.winStarts     : start indices of fixed winSec windows
%   D.winEnds       : end indices of fixed winSec windows
%   D.thr           : numeric threshold on activity
%   D.activity      : activity signal used for detection (x.^4 smoothed)
%   D.mask          : logical mask after cleanup (1=burst, 0=silence)

    if nargin < 3 || isempty(opts), opts = struct(); end
    opts = fillDefaults(opts);

    X = double(X);
    if size(X,2) ~= 2
        error('X must be [N x 2].');
    end
    N = size(X,1);

    Ns = round(opts.winSec * Fs);
    half = floor(Ns/2);

    % Compute activity on both channels and select stronger channel
    a1 = activity_x4(X(:,1), Fs, opts.smoothMs);
    a2 = activity_x4(X(:,2), Fs, opts.smoothMs);

    if max(a2) > max(a1)
        chUsed = 2; activity = a2;
    else
        chUsed = 1; activity = a1;
    end

    thr = opts.thrFrac * max(activity);
    mask = activity >= thr;

    % Cleanup (debounce)
    minOn  = round(opts.minOnMs  * 1e-3 * Fs);
    minOff = round(opts.minOffMs * 1e-3 * Fs);
    mask = removeShortRuns(mask, minOn,  true);
    mask = removeShortRuns(mask, minOff, false);

    % Burst boundaries
    [starts, ends] = runs(mask);

    % If nothing found, return empty arrays but keep fields
    if isempty(starts)
        D = makeOutput(chUsed, [], [], [], [], [], thr, activity, mask);
        return;
    end

    % Keep up to maxReps longest bursts
    lens = ends - starts + 1;
    [~, order] = sort(lens, 'descend');
    order = order(1:min(opts.maxReps, numel(order)));

    starts = starts(order);
    ends   = ends(order);

    % Sort selected bursts by time (nice for plotting/reporting)
    [starts, idx] = sort(starts);
    ends = ends(idx);

    % Midpoints and fixed analysis windows
    midpoints = round((starts + ends)/2);
    winStarts = midpoints - half;
    winEnds   = winStarts + Ns - 1;

    % Clamp windows to valid range
    for i = 1:numel(winStarts)
        if winStarts(i) < 1
            winStarts(i) = 1;
            winEnds(i) = min(N, Ns);
        end
        if winEnds(i) > N
            winEnds(i) = N;
            winStarts(i) = max(1, N - Ns + 1);
        end
    end

    D = makeOutput(chUsed, starts, ends, midpoints, winStarts, winEnds, thr, activity, mask);
end

%% ---------------- Local helper functions ----------------
function opts = fillDefaults(opts)
    if ~isfield(opts,'thrFrac'),  opts.thrFrac = 0.02; end
    if ~isfield(opts,'smoothMs'), opts.smoothMs = 10; end
    if ~isfield(opts,'minOnMs'),  opts.minOnMs = 30; end
    if ~isfield(opts,'minOffMs'), opts.minOffMs = 20; end
    if ~isfield(opts,'winSec'),   opts.winSec = 0.5; end
    if ~isfield(opts,'maxReps'),  opts.maxReps = 4; end
end

function a = activity_x4(x, Fs, smoothMs)
    x = double(x(:));
    x = x - mean(x);
    a = x.^4;
    win = max(1, round((smoothMs/1000)*Fs));
    a = movmean(a, win);
end

function [starts, ends] = runs(logicalVec)
    logicalVec = logical(logicalVec(:));
    d = diff([false; logicalVec; false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
end

function y = removeShortRuns(x, minLen, onValue)
% If onValue==true: remove short ON runs (set them to 0)
% If onValue==false: fill short OFF gaps (set them to 1)
    x = logical(x(:));
    y = x;

    if minLen <= 1
        return;
    end

    if onValue
        [s,e] = runs(x);
        lens = e - s + 1;
        kill = lens < minLen;
        for i = find(kill).'
            y(s(i):e(i)) = false;
        end
    else
        [s,e] = runs(~x);
        lens = e - s + 1;
        fill = lens < minLen;
        for i = find(fill).'
            y(s(i):e(i)) = true;
        end
    end
end

function D = makeOutput(chUsed, starts, ends, midpoints, winStarts, winEnds, thr, activity, mask)
    D = struct();
    D.chUsed    = chUsed;
    D.starts    = starts(:);
    D.ends      = ends(:);
    D.midpoints = midpoints(:);
    D.winStarts = winStarts(:);
    D.winEnds   = winEnds(:);
    D.thr       = thr;
    D.activity  = activity(:);
    D.mask      = mask(:);
end
