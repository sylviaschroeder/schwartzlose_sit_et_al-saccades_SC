function  [shifts, shift_onset_T, shift_lims] = lin_shift(time, onset_T, nShifts, shift_type)

if nargin <4
    shift_type = 'linear';
end

nTs = numel(time);

switch shift_type
    case 'linear'

        shift_lims = round(prctile(onset_T(:,1), [20 80]));
        shift_ratio = shift_lims(1)/(nTs - shift_lims(2));
        min_win = 0.5*nShifts*0.1/mean(diff(time)); % this is hardcoded. Is the minimum required average distance betwen shifts. Now 100ms.
        left_win = max(shift_ratio*min_win, shift_lims(1));
        right_win = max(nTs - shift_lims(2), min_win/shift_ratio);

        shift_neg = floor(rand(floor(0.5*nShifts*left_win/right_win), 1)*(-1)*left_win) +1; % make shift proportional
        shift_pos = floor(rand(floor(0.5*nShifts*right_win/left_win), 1)*(right_win)) - 1;
        shifts = [0; shift_neg; shift_pos];
        shifts = shifts(1:nShifts);

        shift_onset_T = onset_T;
        shift_onset_T(onset_T< left_win | onset_T >= nTs-right_win) =[];
        shift_onset_T = shift_onset_T + shifts'; % nSaccade * (nShifts+1)
        shift_lims = [left_win, nTs-right_win];

    case 'circ'

        max_shift = floor(nTs/2);
        shift_neg = round(rand(round(nShifts/2),1)*max_shift);
        shift_pos = round(rand(round(nShifts/2),1)*(-1)*max_shift);
        shifts = [0; shift_neg; shift_pos];
        shifts = shifts(1:nShifts);

        shift_onset_T = onset_T;
        shift_onset_T = shift_onset_T + shifts'; % nSaccade * (nShifts+1)
        shift_onset_T(shift_onset_T>nTs) = shift_onset_T(shift_onset_T>nTs) -nTs;
        shift_onset_T(shift_onset_T<=0) = shift_onset_T(shift_onset_T<=0) + nTs;

        shift_lims = [0, nTs];

end

end