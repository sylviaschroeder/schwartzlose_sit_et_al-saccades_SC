function [x, fitResp] = fit2dGaussRF_wShift(rf, makeplots)
% resp = your measure of the responsiveness of the neuron at each point in
% space
% x = [amplitude, xCenter, xStdev, yCenter, yStdev, rotation, xShift2, xShift3, ...]
    rf_all = mean(rf,3);

    xcoords = 1:size(rf,2);
    ycoords = 1:size(rf,1);

    [x, y] = meshgrid(xcoords, ycoords);
    xdata = zeros(size(x,1), size(x,2), 2);
    xdata(:,:,1) = x;
    xdata(:,:,2) = y;
    
    % x0 is the initial guess
    [maxY, maxX] = find(rf_all == max(rf_all(:)), 1);
    x0 = [max(rf_all(:)), maxX, mean(diff(xcoords))*5, maxY, ...
        mean(diff(ycoords))*5, 0, zeros(1, size(rf,3)-1)];        

    lb = [0, min(xcoords), 0, min(ycoords), 0, -pi/4, ...
        -ones(1, size(rf,3)-1) .* max(xcoords)];
    mx_std = max(range(xcoords(:)), range(ycoords(:))) / 4;
    ub = [2 * max(rf_all(:)), max(xcoords), mx_std, ...
        max(ycoords), mx_std, pi/4, ...
        ones(1, size(rf,3)-1) .* max(xcoords)];
    options = optimoptions('lsqcurvefit', 'Display', 'off');
    x = lsqcurvefit(@whiteNoise.D2GaussFunctionRot_wShift, x0, xdata, ...
        rf, lb, ub, options);

    fitResp = whiteNoise.D2GaussFunctionRot_wShift(x, xdata);
    
    if makeplots
        numShifts = size(rf,3);
        figure
        for sh = 1:numShifts
            subplot(3,numShifts,sh);
            imagesc(xcoords, ycoords, rf(:,:,sh));
            colormap hot;

            subplot(3,numShifts,numShifts+sh);
            imagesc(xcoords, ycoords, fitResp(:,:,sh));

            subplot(3,numShifts,2*numShifts+sh);
            contour(fitResp(:,:,sh),1,'Color', 'k', 'LineWidth', 2.0);
            set(gca, 'YDir', 'reverse');
        end
    end
end