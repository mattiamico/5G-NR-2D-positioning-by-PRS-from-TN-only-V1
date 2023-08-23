function rstd = getRSTDValues(toa, sr)
%   RSTD = getRSTDValues(TOA,SR) computes RSTD values, given the time of
%   arrivals TOA and the sampling rate SR.

    % Compute number of samples delay between arrivals from different gNBs
    rstd = zeros(length(toa));
    for jj = 1:length(toa)
        for ii = 1:length(toa)
            rstd(ii,jj) = toa(ii) - toa(jj);
        end
    end

    % Get RSTD values in time
    rstd = rstd./sr; % [s]
end

