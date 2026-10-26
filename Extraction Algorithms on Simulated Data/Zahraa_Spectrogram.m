%% Radar spectrogram
% Zahraa Sattar | EEE4022S
%
% Calculate the windowed Fourier transform of the radar signal.
% Return the complex spectrum with its frequency and time axes.

function [S, freq_axis, time_axis] = Zahraa_Spectrogram(x, win, noverlap, nfft, fs)

    x   = x(:).';
    win = win(:).';
    L = length(win);        % frame length = window length
    N = length(x);
    hop = L - noverlap;     % step between successive frame starts

    M = floor((N - L) / hop) + 1;    % number of frames
    S = zeros(nfft, M);

    frame_starts = zeros(1, M);      % keep track of where each frame starts (for time axis later)

    for m = 1:M
        start_idx = (m - 1) * hop + 1;
        stop_idx  = start_idx + L - 1;
        frame_starts(m) = start_idx;
        segment = x(start_idx:stop_idx) .* win;     % grab the chunk and window it
        S(:, m) = fft(segment, nfft).';
    end

    time_axis = (frame_starts - 1 + L/2) / fs;  % labelling the time_axis as the center of the window

    if isreal(x)
        freq_axis = (0:nfft-1) * (fs / nfft);
    else
        S = fftshift(S, 1);
        freq_axis = (-nfft/2 : nfft/2 - 1) * (fs / nfft);
    end
end
