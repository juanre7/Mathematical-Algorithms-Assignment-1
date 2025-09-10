Here’s a **cleaned, corrected, and tightened** version of your README. I fixed parameter inconsistencies, corrected the sampling period, made quantization well-defined (uniform mid-rise over $[-1,1]$), and added small clarity tweaks.

# From Analog to Digital: Signal Simulation

![Header Image](https://github.com/juanre7/Mathematical-Algorithms-Assignment-1/blob/main/header.png?raw=true)

## Introduction

This MATLAB walkthrough shows the pipeline from a continuous-time signal to a digital bitstream:
1) generate an “analog” reference signal (dense grid),  
2) sample it,  
3) uniformly quantize the samples,  
4) encode to binary, and  
5) concatenate into a bitstream.

Concepts used include the Nyquist criterion (sample ≥ 2× the highest frequency) and aliasing (distortion when sampling below that). The code is self-contained so you can vary parameters and immediately see the effect.

### Learning Objectives
- Understand the full path from analog signals to digital representation.
- Simulate each step in MATLAB: generation, sampling, quantization, encoding.
- Connect these ideas to practical sensor/IoT acquisition.

---

## Step 0: Simulation Parameters (centralized)

```matlab
% Core parameters
f0   = 1000;          % Signal frequency in Hz
Tend = 0.01;          % Observation window = 10 ms (~10 periods at 1 kHz)

% Dense "analog" grid to approximate continuity
t        = 0:1e-6:Tend;          % 1 µs step
x_analog = sin(2*pi*f0*t);
````

---

## Step 1: Generate the Analog Signal

```matlab
figure;
plot(t, x_analog, 'LineWidth', 1.5);
grid on;
title('Analog Reference Signal (Sine, 1 kHz)');
xlabel('Time (s)'); ylabel('Amplitude');
```

This is a pure 1 kHz sine, a simple stand-in for vibration/audio sensor data.

---

## Step 2: Sampling

```matlab
Fs = 5000;                  % Sampling frequency (5 kHz > 2*f0 to avoid aliasing)
Ts = 1/Fs;                  % Sampling period (correct)
n  = 0:Ts:Tend;             % Sample instants
x_samp = sin(2*pi*f0*n);    % Sampled values

figure;
stem(n, x_samp, 'filled'); hold on;
plot(t, x_analog, 'LineWidth', 1.0);  % overlay for visual comparison
grid on;
title(sprintf('Sampling at Fs = %.0f Hz', Fs));
xlabel('Time (s)'); ylabel('Amplitude');
legend('Samples','Analog (overlay)','Location','best');
```

> Tip: set `Fs` below, at, and above `2*f0` to observe aliasing.

---

## Step 3: Quantization (Uniform, Mid-Rise)

We quantize to a fixed, known range $[-1,1]$ that covers the sine.
(Mid-rise: reconstruction levels sit at the center of each bin; for even `L`, zero is *between* two levels.)

```matlab
bits = 4;                    % Number of bits per sample
L    = 2^bits;               % Number of quantization levels

xmin = -1; xmax = 1;         % Quantizer range (covers the sine)
Delta = (xmax - xmin)/L;     % Step size

% Clip to range to avoid out-of-range indices (e.g., numeric edge cases)
x_clip = min(max(x_samp, xmin), xmax - eps);

% Mid-rise index and reconstruction value
idx = floor((x_clip - xmin)/Delta);     % indices: 0..L-1
x_q = xmin + (idx + 0.5)*Delta;         % reconstructed amplitudes

% Metrics
e   = x_q - x_samp;                      % quantization error
SQNR_dB = 10*log10(var(x_samp)/var(e));  % signal-to-quantization-noise ratio
MSE     = mean(e.^2);

figure;
stairs(n, x_q, 'LineWidth', 1.25); hold on;
stem(n, x_samp, 'filled');
grid on;
title(sprintf('Uniform Quantization: %d levels (%d bits) | SQNR = %.2f dB | MSE = %.3e', L, bits, SQNR_dB, MSE));
xlabel('Time (s)'); ylabel('Amplitude');
legend('Quantized (stairs)','Original samples','Location','best');
```

> Rule of thumb for a full-scale sine: `SQNR_dB ≈ 6.02·bits + 1.76`.

---

## Step 4: Binary Encoding

Convert indices `idx ∈ {0,…,L−1}` to fixed-length binary words.

```matlab
binary_words = dec2bin(idx, bits);   % char array (N_samples × bits)

disp('--- First 10 encoded samples (binary words) ---');
disp(binary_words(1:min(10, size(binary_words,1)), :));
```

---

## Step 5: Digital Bitstream

Concatenate all codewords into a single stream (string of `'0'`/`'1'`).
Optionally, also produce a numeric vector of 0/1.

```matlab
% As a single character string
bitstream_char = reshape(binary_words.', 1, []);   % column-wise to row
disp('--- First 40 bits of the stream (char) ---');
disp(bitstream_char(1:min(40, numel(bitstream_char))));

% As a numeric vector (uint8 of 0/1)
bitstream_num = uint8(bitstream_char.' - '0');     % column vector of 0/1
```

---

## Summary

```matlab
fprintf('\nSimulation complete!\n');
fprintf('Analog -> Sampling -> Quantization -> Binary Encoding -> Bitstream\n');
fprintf('Bits per sample: %d\n', bits);
fprintf('Samples: %d\n', numel(x_samp));
fprintf('Total bits: %d\n', numel(bitstream_num));
fprintf('SQNR: %.2f dB | MSE: %.3e\n', SQNR_dB, MSE);
```

---

## Homework Exercises

### Exercise 1 — Sampling vs Nyquist

![Header Image](https://raw.githubusercontent.com/juanre7/Mathematical-Algorithms-Assignment-1/refs/heads/main/Nyquist.bmp)


```matlab
f0 = 1000; Tend = 0.01; t = 0:1e-6:Tend; x_analog = sin(2*pi*f0*t);
Fs_list = [1500, 2000, 5000];  % <Nyquist, =Nyquist, >Nyquist

figure('Name','Sampling Comparison','NumberTitle','off');
for k = 1:numel(Fs_list)
    Fs = Fs_list(k); Ts = 1/Fs; n = 0:Ts:Tend;
    x_s = sin(2*pi*f0*n);

    subplot(numel(Fs_list),1,k);
    plot(t, x_analog, 'LineWidth', 1.0); hold on;
    stem(n, x_s, 'filled'); grid on;
    title(sprintf('Fs = %.0f Hz', Fs));
    xlabel('Time (s)'); ylabel('Amplitude');

    if Fs < 2*f0
        text(0.001, -0.85, 'Aliasing expected (< 2f0)', 'FontWeight', 'bold');
    elseif Fs == 2*f0
        text(0.001, -0.85, 'At Nyquist: phase/noise sensitive', 'FontWeight', 'bold');
    else
        text(0.001, -0.85, 'Above Nyquist: safe sampling', 'FontWeight', 'bold');
    end
end
```

### Exercise 2 — Quantization Levels (8, 16, 64)

![Exercise 2](https://raw.githubusercontent.com/juanre7/Mathematical-Algorithms-Assignment-1/refs/heads/main/quantization.bmp)


```matlab
% Safe sampling to isolate quantization effects
f0=1000; Tend=0.01; Fs=5*f0; Ts=1/Fs; n=0:Ts:Tend;
x_samp = sin(2*pi*f0*n);
xmin=-1; xmax=1;

L_list = [8, 16, 64];

figure('Name','Quantization Levels','NumberTitle','off');
for k = 1:numel(L_list)
    L=L_list(k); bits=log2(L); Delta=(xmax-xmin)/L;
    x_clip = min(max(x_samp, xmin), xmax - eps);
    idx = floor((x_clip - xmin)/Delta);
    x_q = xmin + (idx + 0.5)*Delta;

    e = x_q - x_samp;
    SQNR = 10*log10(var(x_samp)/var(e));

    subplot(numel(L_list),1,k);
    stairs(n, x_q, 'LineWidth', 1.2); hold on; stem(n, x_samp,'filled');
    grid on; title(sprintf('L=%d (%.0f bits), \\Delta=%.3f, SQNR=%.2f dB', L, bits, Delta, SQNR));
    xlabel('Time (s)'); ylabel('Amplitude');
    legend('Quantized','Original','Location','best');
end
```

---

## Notes

* **Nyquist criterion:** to avoid ambiguity, sample at or above twice the highest frequency present.
* **Aliasing:** when sampled below Nyquist, higher-frequency content “folds” to lower frequencies.
* **Uniform quantization (mid-rise):** step size $\Delta=(x_{\max}-x_{\min})/L$; error is bounded by $\pm\Delta/2$.
* **SQNR (Signal-to-Quantization-Noise Ratio):** for a full-scale sine, roughly $6.02\,\text{dB/bit} + 1.76$.

---

*Last updated: 2025-09-10*
*Author: \Juan Rodriguez Esteban*

```
```
