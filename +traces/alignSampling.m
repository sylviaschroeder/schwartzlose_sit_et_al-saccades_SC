function [newTraces, newTime] = alignSampling(traces, time, delayIDs, ...
    delays)

timeBin = median(diff(time));

if length(delays) <= 1
    newTraces = traces;
    newTime = time;
    return
end

delta_t = median(diff(delays));
upsample = round(timeBin / delta_t);
timeBin = timeBin / upsample;
newTime = reshape((time + (0:upsample-1) * timeBin)', [], 1);
newTraces = NaN(length(newTime), size(traces,2));
for d = 1:length(delays)
    indUnits = find(delayIDs == d);
    for n = indUnits'
        if sum(isnan(traces(:,n)))/size(traces,1) > 0.8
            continue
        end
        nanInd1 = isnan(traces(:,n));
        newTraces(:,n) = interp1(time(~nanInd1) + delays(d), ...
            traces(~nanInd1,n), newTime, 'pchip');
        nanInd2 = reshape(repmat(nanInd1, 1, upsample)', [], 1);
        newTraces(nanInd2,n) = NaN;
    end
end