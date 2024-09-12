function [stimMatrix, x, y, z] = ...
    getStimMatrix(times, xyPos, diameter, isWhite)

% determine tested circle positions and sizes
x = unique(xyPos(:,1));
y = unique(xyPos(:,2));
z = unique(diameter);

% check whether all combinations of positions and sizes were presented; if
% not, set stimMatrix for those combinations to NaN
allCombis = [ ...
    reshape(repmat(x(:), 1, length(y), length(z)), [], 1), ...
    reshape(repmat(y(:)', length(x), 1, length(z)), [], 1), ...
    reshape(repmat(permute(z(:),[2 3 1]), ...
    length(x), length(y), 1), [], 1)];
notTested = setdiff(allCombis, unique([xyPos, diameter], 'rows'), 'rows');

% for each time point (1st dim), make grid of all tested positions (2nd + 
% 3rd dim), and all tested sizes (4th dim); 1: white circle, -1: black
% circle
stimMatrix = zeros(length(times), length(x), length(y), length(z));
% row indices for each time
r = xyPos(:,1)' == x;
indR = repmat((1:length(x))', 1, length(times));
r = indR(r);
% column indices for each time
c = xyPos(:,2)' == y;
indC = repmat((1:length(y))', 1, length(times));
c = indC(c);
% diameter indices for each time
d = diameter' == z;
indD = repmat((1:length(z))', 1, length(times));
d = indD(d);
% get linear indices of circle positions and sizes for each time
ind = sub2ind(size(stimMatrix), (1:length(times))', r, c, d);
stimMatrix(ind(isWhite == 1)) = 1;
stimMatrix(ind(isWhite ~= 1)) = -1;

for j = 1:size(notTested,1)
    stimMatrix(:,notTested(1),notTested(2),notTested(3)) = NaN;
end