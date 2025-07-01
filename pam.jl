using TyPlot, FFTW

# 参数设置
signal_freq = 1000       # 信号频率 (Hz)
fs = 8000                # 采样频率 (Hz)
duty_cycle = 0.8         # 脉冲占空比
oversample_factor = 200  # 过采样因子
num_cycles = 10          # 模拟周期数

# 时间轴
T_signal = 1 / signal_freq
total_time = num_cycles * T_signal
t = range(0, total_time, length=round(Int, total_time * fs * oversample_factor))

# 1. 生成正弦波信号
function generate_sine_wave(t, f=signal_freq)
    return sin.(2π * f * t)
end

# 生成方波信号
function generate_square_wave(t, f=signal_freq)
    return sign.(sin.(2π * f * t))  # sign函数将正弦波转换为±1的方波
end

# 生成锯齿波信号
function generate_sawtooth_wave(t, f=signal_freq)
    return 2 * mod.(f * t, 1.0) .- 1  # 将时间映射到[-1,1]的锯齿波
end


# 正弦信号生成函数
generate_sine_wave(t, f=signal_freq) = sin.(2π * f * t)

# 方波信号生成函数
#generate_sine_wave(t, f=signal_freq) = sin.(2π * f * t)

# 锯齿信号生成函数
#generate_sine_wave(t, f=signal_freq) = sin.(2π * f * t)

# 曲顶采样（自然采样）
function generate_natural_sampling_pulses(t, fs, signal, duty_cycle)
    pulse_period = 1 / fs
    pulse_width = duty_cycle * pulse_period
    pulses = zeros(length(t))
    sample_times = 0:pulse_period:maximum(t)

    for i in eachindex(t)
        nearest_sample_time = sample_times[argmin(abs.(t[i] .- sample_times))]
        if abs(t[i] - nearest_sample_time) <= pulse_width / 2
            pulses[i] = signal[i]
        end
    end
    return pulses
end

# 汉宁窗（用于频谱可视化）
function hanning(N)
    n = 0:N-1
    return 0.5 .* (1 .- cos.(2π * n / (N-1)))
end

# 频谱计算
function compute_spectrum(signal, fs)
    N = length(signal)
    window = hanning(N)
    windowed_signal = signal .* window
    fft_result = fft(windowed_signal)
    fft_magnitude = abs.(fft_result) / N * 2
    freq = (0:N-1) .* fs / N
    half_idx = Int(floor(N/2))
    return freq[1:half_idx], fft_magnitude[1:half_idx]
end

# 频域低通滤波重建
function frequency_domain_lowpass_filter_reconstruction(signal::Vector{Float64}, fs::Float64, cutoff_freq::Float64)
    N = length(signal)
    fft_result = fft(signal)
    freq = collect(0:N-1) .* fs / N
    freq_shifted = fftshift(freq .- fs/2)
    H = zeros(ComplexF64, N)
    for i in 1:N
        if abs(freq_shifted[i]) <= cutoff_freq
            H[i] = 1.0
        end
    end
    fft_shifted = fftshift(fft_result)
    filtered_shifted = fft_shifted .* H
    filtered_fft = ifftshift(filtered_shifted)
    reconstructed_signal = real(ifft(filtered_fft))
    return reconstructed_signal
end

# 主函数
function plot_all()
    input_signal = generate_sine_wave(t)
    pulses = generate_natural_sampling_pulses(t, fs, input_signal, duty_cycle)
    freq, spectrum = compute_spectrum(pulses, fs)
    reconstructed_signal = frequency_domain_lowpass_filter_reconstruction(pulses, float(fs), 5000.0)

    display_range = findall(t .< 2 * T_signal)

    # 原始信号、PAM、频谱图
    figure(figsize=(12, 10))

    subplot(3, 1, 1)
    plot(t[display_range], input_signal[display_range], "b-", linewidth=1.5)
    title("原始正弦信号 ($(signal_freq)Hz)")
    ylabel("幅度")
    grid(true)

    subplot(3, 1, 2)
    plot(t[display_range], pulses[display_range], "r-", linewidth=1.0)
    title("曲顶采样信号 (fs=$(fs)Hz，占空比=$(round(duty_cycle * 100, digits=1))%)")
    ylabel("幅度")
    grid(true)

    subplot(3, 1, 3)
    plot(freq, spectrum, "g-", linewidth=1.5)
    xlim(0, 5000)
    title("频谱 (抽样信号)")
    xlabel("频率 (Hz)")
    ylabel("幅度")
    grid(true)

    # 重建图
    figure(2)
    plot(t[display_range], reconstructed_signal[display_range], "m-", linewidth=2)
    title("频域低通滤波重建信号 (cutoff=3000Hz)")
    xlabel("时间 (s)")
    ylabel("幅度")
    grid(true)
end

# 执行
plot_all()
