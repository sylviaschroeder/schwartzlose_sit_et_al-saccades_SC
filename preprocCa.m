function newData = preprocCa(caData, opts)

% for correcting baseline drifts of calcium traces at start of experiments
driftWin = 20; % in s, window to test whether baseline is higher than normal
driftThresh = 1.5; % in std, threshold for drift
correctWin = 150; % in s, window to fit exponential


    %% Prepare calcium traces
        % interpolate calcium traces to align all to same time
        t_ind = caData.time;
        traces = caData.traces;
        t_ca = caData.time;
        timeBin = median(diff(t_ca));
        if length(caData.delays) > 1
            delta_t = median(diff(caData.delays));
            upsample = round(timeBin / delta_t);
            timeBin = timeBin / upsample;
            t_up = reshape((t_ca + (0:upsample-1) * timeBin)', [], 1);
            tr_up = NaN(length(t_up), size(traces,2));
            for d = 1:length(caData.delays)
                indUnits = find(caData.planes == d);
                for n = indUnits'
                    if all(isnan(traces(:,n)))
                        continue
                    end
                    nanInd1 = isnan(traces(:,n));
                    tr_up(:,n) = interp1(t_ca(~nanInd1) + caData.delays(d), ...
                        traces(~nanInd1,n), t_up, 'pchip');
                    nanInd2 = reshape(repmat(nanInd1, 1, upsample)', [], 1);
                    tr_up(nanInd2,n) = NaN;
                end
            end
            t_ca = t_up;
            traces = tr_up;
        end

        % remove strong baseline decay at start of experiment in cells that
        % show it
        indUnits = find(mean(traces(1:round(driftWin / timeBin),:), 1, 'omitnan') > ...
            mean(traces, 1, 'omitnan') + driftThresh .* std(traces,0,1, 'omitnan'));
        ind = round(correctWin / timeBin);
        for iUnit = 1:length(indUnits)
            y = traces(1:ind, indUnits(iUnit));
            y = fillmissing(y, 'linear');
            % fit double exponential to start of trace
            f = fit((1:length(y))', y, ...
                @(a,b,c,d,e,x) a + b .* exp(-x ./ c) + d .* exp(-x ./ e), ...
                'Lower', [0 0 0 0 0], ...
                'Upper', [max(y) max(y) 500 max(y) 500], ...
                'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
            % remove fit
            traces(:, indUnits(iUnit)) = traces(:, indUnits(iUnit)) - ...
                f(1 : size(traces,1)) + f.a;
        end

        newData.traces = traces;
        newData.time = t_ca;
        newData.ids = caData.ids;
        newData.planes = caData.planes;
end