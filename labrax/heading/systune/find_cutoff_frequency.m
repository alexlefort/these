function wc = find_cutoff_frequency(sys)
    [mag,~,wout] = bode(sys);                                   % Get Plot Data
    mag = squeeze(mag);                                             % Reduce (1x1xN) Matrix To (1xN)
    magr2 = (mag).^2;                                      % Calculate Power Of Ratio Of ‘mag/max(mag)’
    wc = interp1(magr2, wout, 0.5, 'spline');