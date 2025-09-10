% Learning Objectives
% 
% Understand the complete path from analog signals to digital representation.
% Simulate each step in MATLAB: signal generation, sampling, quantization, and encoding.
% Relate the concepts to IoT sensor data acquisition.

%% 1. Generate Analog Signal (conceptual) - generate a simple continuous-time signal (e.g., sine wave).

t = 0:0.0001:0.01;       % very fine step (continuous-like time)

f = 1000;                 % frequency = 100 Hz

x_analog = sin(2*pi*f*t);

 

figure;

plot(t, x_analog, 'LineWidth', 1.5);

title('Analog Signal (Sine Wave)');

xlabel('Time (s)'); ylabel('Amplitude');

 


%% 2. Sampling - Sample the analog signal at a chosen sampling frequency.


Fs = 1000;               % Sampling frequency = 1 kHz

Ts = 0.1/Fs;               % Sampling period

n = 0:Ts:0.01;           % Discrete sample points

x_sampled = sin(2*pi*f*n);

 

figure;

stem(n, x_sampled, 'filled');

title('Sampled Signal');

xlabel('Time (s)'); ylabel('Amplitude');

 

 % Discussion: Compare analog vs sampled wave. Relate to Nyquist theorem.



%% 3. Quantization - Simulate quantization by rounding amplitudes to finite levels.


bits = 4;                           % Number of bits

levels = 2^bits;                    % Quantization levels

x_min = min(x_sampled);

x_max = max(x_sampled);

q_step = (x_max - x_min)/levels;    % Step size

 

x_index = round((x_sampled - x_min)/q_step); % Map samples to indices

x_quantized = x_index*q_step + x_min;        % Map back to amplitude

 

figure;

stem(n, x_quantized, 'filled');

title(['Quantized Signal (' num2str(bits) '-bit)']);

xlabel('Time (s)'); ylabel('Amplitude');

 

 % Discussion: More levels → higher accuracy, less compression.



%% 4. Encoding (Binary Representation) - Convert quantized values into binary codes.


binary_codes = dec2bin(x_index, bits); % Convert indices to binary words

 

disp('--- First 10 encoded samples ---');

disp(binary_codes(1:10,:));

 

 

% Now see how each quantized sample becomes a digital word.

 



%% 5. Digital Bitstream - Show the digital stream (bit sequence).

bitstream = reshape(binary_codes.',1,[]); % Concatenate into one string

disp('--- First 40 bits of the stream ---');

disp(bitstream(1:40));

 

 

%% Summary

fprintf('\nSimulation complete!\n');

fprintf('Analog -> Sampling -> Quantization -> Binary Encoding -> Digital Stream\n');

fprintf('Bits per sample: %d\n', bits);

fprintf('Total samples: %d\n', length(x_sampled));

fprintf('Total bits: %d\n', length(bitstream));



%% Homework
% 
% Modify the MATLAB script:
%% Try different sampling frequencies (below Nyquist, at Nyquist, and above Nyquist).

% Reference "analog" signal (very densely sampled)
f0     = 1000;           % Hz (signal frequency)
Tend   = 0.01;           % 10 ms observation (~10 periods at 1 kHz)
t      = 0:1e-6:Tend;    % very fine step (1 µs)
x_anal = sin(2*pi*f0*t);

% Sampling frequency sets
Fs_cases = [1.5*f0, 2*f0, 5*f0];    % <Nyquist, =Nyquist, >Nyquist
labels   = {"Below Nyquist (Fs = 1.5·f0)", ...
            "At Nyquist (Fs = 2·f0)", ...
            "Above Nyquist (Fs = 5·f0)"};   

Nyq = 2*f0;  % Nyquist frequency (Hz)

fprintf('Signal frequency f0 = %.0f Hz | Nyquist = %.0f Hz\n', f0, Nyq);
fprintf('Fs cases: %s\n', mat2str(Fs_cases));

figure('Name','Sampling comparison','NumberTitle','off');

for k = 1:numel(Fs_cases)
    Fs = Fs_cases(k);
    Ts = 1/Fs;
    n  = 0:Ts:Tend;
    xS = sin(2*pi*f0*n);

    % Informative calculation of the expected alias (in Hz) for Fs < 2*f0
    % f_alias = |f0 - m*Fs| with m = round(f0/Fs), constrained to [0, Fs/2]
    m = round(f0/Fs);
    f_alias = abs(f0 - m*Fs);
    if f_alias > Fs/2
        f_alias = Fs - f_alias;  % fold into the half-band [0, Fs/2]
    end

    subplot(3,1,k);
    plot(t, x_anal, 'LineWidth', 1.25); hold on;
    stem(n, xS, 'filled'); % samples
    grid on;
    title(sprintf('%s  |  Fs = %.0f Hz', labels{k}, Fs));
    xlabel('Time (s)'); ylabel('Amplitude');
    legend('Continuous signal (fine)', 'Samples', 'Location', 'best');

    if Fs < Nyq
        txt = sprintf('ALIASING: f_{alias} ≈ %.1f Hz (expected)', f_alias);
        text(0.001, -0.8, txt, 'FontWeight','bold');
    elseif Fs == Nyq
        text(0.001, -0.8, 'At Nyquist: ambiguous phase/amplitude detection', 'FontWeight','bold');
    else
        text(0.001, -0.8, 'Above Nyquist: safe sampling', 'FontWeight','bold');
    end

    % Console message
    if Fs < Nyq
        fprintf('Fs = %.0f Hz (< Nyquist): expected alias ≈ %.1f Hz\n', Fs, f_alias);
    elseif Fs == Nyq
        fprintf('Fs = %.0f Hz (= Nyquist): edge case; sensitive to phase/noise\n', Fs);
    else
        fprintf('Fs = %.0f Hz (> Nyquist): proper sampling\n', Fs);
    end
end

%% Use different quantization levels (8, 16, 64). Plot and compare the effect on signal quality.
f0     = 1000;             % Hz
Tend   = 0.01;             % 10 ms
t      = 0:1e-6:Tend;      % dense "analog" grid
x_anal = sin(2*pi*f0*t);

Fs  = 5*f0;                % safe sampling to isolate quantization
Ts  = 1/Fs;
n   = 0:Ts:Tend;
xS  = sin(2*pi*f0*n);

L_set = [8, 16, 64];
xmin = -1; xmax = 1;

fprintf('Fs = %.0f Hz (Nyquist = %.0f Hz)\n', Fs, 2*f0);
fprintf('L   Bits  Delta      SQNR(dB)  MSE\n');

figure('Name','Quantization (overview + zoom + error)','NumberTitle','off');

for k = 1:numel(L_set)
    L = L_set(k);
    B = log2(L);
    Delta = (xmax - xmin)/L;

    % Mid-rise uniform quantizer
    x_clip = min(max(xS, xmin), xmax - eps);
    idx = floor((x_clip - xmin)/Delta);
    xQ  = xmin + (idx + 0.5)*Delta;

    % Metrics
    e    = xQ - xS;
    SQNR = 10*log10(var(xS)/var(e));
    MSE  = mean(e.^2);

    fprintf('%-3d %-4.0f  %-9.5f  %8.2f  % .3e\n', L, B, Delta, SQNR, MSE);

    % -------- Overview: full 10 ms with quantization levels
    subplot(3,3,3*(k-1)+1);
    plot(t, x_anal, 'LineWidth', 1.0); hold on;
    stairs(n, xQ, 'LineWidth', 1.25);
    ylv = linspace(xmin+Delta/2, xmax-Delta/2, L); % level centers
    yline(ylv, ':');                                % show levels
    grid on; axis([0 Tend xmin xmax]);
    title(sprintf('L=%d (%.0f bits) | Δ=%.3f | SQNR=%.2f dB', L, B, Delta, SQNR));
    xlabel('Time (s)'); ylabel('Amplitude');
    legend('Continuous','Quantized','Levels','Location','southoutside');

    % -------- Zoomed view: first 2 ms to make steps obvious
    subplot(3,3,3*(k-1)+2);
    plot(t, x_anal, 'LineWidth', 1.0); hold on;
    stairs(n, xQ, 'LineWidth', 1.25);
    yline(ylv, ':');
    grid on; axis([0 2e-3 xmin xmax]);
    title('Zoom: first 2 ms');
    xlabel('Time (s)'); ylabel('Amplitude');

    % -------- Quantization error vs time
    subplot(3,3,3*(k-1)+3);
    stem(n, e, 'filled'); grid on;
    axis([0 2e-3 -Delta Delta]); % scale to ±Δ
    title(sprintf('Error e[n] (±Δ/2 ≈ ±%.3f)', Delta/2));
    xlabel('Time (s)'); ylabel('Error');
end

% Theoretical SQNR for a full-scale sine (for reference)
B_theory = log2(L_set);
SQNR_theory = 6.02*B_theory + 1.76;
disp(table(L_set.', B_theory.', SQNR_theory.', 'VariableNames', ...
    {'L','Bits','SQNR_theory_dB'}));
