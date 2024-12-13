function [rf_visDeg, x1, y1] = interpolateRFtoVisualDegrees(rf, stimEdges)

gridW = diff(stimEdges(1:2)) / size(rf,2);
gridH = -diff(stimEdges(3:4)) / size(rf,1);

% gridpoints of rf in visual degrees
gridX = linspace(stimEdges(1)+0.5*gridW, stimEdges(2)-0.5*gridW, size(rf,2));
gridY = linspace(stimEdges(3)-0.5*gridH, stimEdges(4)+0.5*gridH, size(rf,1));
[x0, y0] = meshgrid(gridX, gridY);

% vectors x1 and y1 specify gridlines with distance of 1 degree
x1 = ceil(gridX(1)) : floor(gridX(end));
y1 = floor(gridY(1)) : -1 : ceil(gridY(end));

rf_visDeg = NaN(size(x1,1), size(x1,2), size(rf,3));
for k = 1:size(rf,3)
    rf_visDeg(:,:,k) = interp2(x0, y0, rf(:,:,k), x1, y1);
end


% % gridpoints of rf
% [x0, y0] = meshgrid((1:size(rf,2))-0.5, (1:size(rf,1))-0.5);
% % vectors x1 and y1 specify gridlines with distance of 1
% % degree (diff(stimPos(...))); values match position of
% % gridlines (edges) in pixels of stimulus row/column
% x1 = linspace(0, size(rf,2), diff(stimEdges(1:2)));
% y1 = linspace(0, size(rf,1), -diff(stimEdges(3:4)));
% 
% % need to delete pixel values outside  given stimulus pixels,
% % so we can use interpolation (rather than extrapolation) when
% % mapping the RF from stimulus pixels to visual degrees
% x2 = x1;
% x2(x1<x0(1) | x1>x0(end)) = [];
% y2 = y1;
% y2(y1<y0(1) | y1>y0(end)) = [];
% [x2, y2] = meshgrid(x2, y2);
% 
% rf_visDeg = NaN(size(x2,1), size(x2,2), size(rf,3));
% for k = 1:size(rf,3)
%     rf_visDeg(:,:,k) = interp2(x0, y0, rf(:,:,k), x2, y2);
% end