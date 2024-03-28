function [x, fitResp] = fit2dGaussRF(rf, makeplots)
% resp = your measure of the responsiveness of the neuron at each point in
% space
% x = [amplitude, xCenter, xStdev, yCenter, yStdev, rotation, baseline]
    xcoords = 1:size(rf,2);
    ycoords = 1:size(rf,1);

    [x, y] = meshgrid(xcoords, ycoords);
    xdata = zeros(size(x,1), size(x,2), 2);
    xdata(:,:,1) = x;
    xdata(:,:,2) = y;
    
    % x0 is the initial guess
    [maxY, maxX] = find(rf == max(rf(:)), 1);
    x0 = [1, maxX, mean(diff(xcoords))*5, maxY, mean(diff(ycoords))*5, 0];
    
    lb = [0, min(xcoords), 0, min(ycoords), 0, -pi/4];
    ub = [2 * max(rf(:)), max(xcoords), (max(xcoords))^2, ...
        max(ycoords), (max(ycoords))^2, pi/4];
    options = optimoptions('lsqcurvefit', 'Display', 'off');
    x = lsqcurvefit(@whiteNoise.D2GaussFunctionRot, x0, xdata, ...
        rf, lb, ub, options);

    fitResp = whiteNoise.D2GaussFunctionRot(x, xdata);
    
    if makeplots
        figure;
        subplot(3,1,1);
        imagesc(xcoords, ycoords, rf);
        colormap hot;
        
        subplot(3,1,2);
        imagesc(xcoords, ycoords, fitResp);
        
        subplot(3,1,3);
        contour(fitResp,1,'Color', 'k', 'LineWidth', 2.0);
        set(gca, 'YDir', 'reverse');
    end
end