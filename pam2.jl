using TyPlot, FFTW

# 参数设置
signal_freq = 1000    # 信号频率 (Hz)
fs = 8000             # 采样频率 (Hz)
duty_cycle = 0.1      # 脉冲占空比
oversample_factor = 200 # 过采样因子
num_cycles = 10       # 分析多个周期以改善频谱分辨率

# 计算时间轴（多个信号周期）
T_signal = 1 / signal_freq  # 信号周期
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

# 2. 真正的曲顶抽样实现
function generate_natural_sampling_pulses(t, fs, signal, duty_cycle)
    pulse_period = 1 / fs  # 采样周期
    pulse_width = duty_cycle * pulse_period
    
    pulses = zeros(length(t))
    
    # 计算所有采样时刻
    sample_times = 0:pulse_period:maximum(t)
    
    for i in eachindex(t)
        # 找到最近的采样时刻
        nearest_sample_time = sample_times[argmin(abs.(t[i] .- sample_times))]
        
        # 计算到最近采样时刻的时间差
        time_diff = abs(t[i] - nearest_sample_time)
        
        # 如果时间点在采样脉冲范围内
        if time_diff <= pulse_width / 2
            # 使用当前时间点的信号值
            pulses[i] = signal[i]
        end
    end
    
    return pulses
end

# 3. 汉宁窗函数
function hanning(N)
    n = 0:N-1
    return 0.5 .* (1 .- cos.(2π * n / (N-1)))
end

# 4. 计算信号的频谱
function compute_spectrum(signal, fs)
    N = length(signal)
    
    # 应用汉宁窗减少频谱泄露
    window = hanning(N)
    windowed_signal = signal .* window
    
    # 计算FFT
    fft_result = fft(windowed_signal)
    fft_magnitude = abs.(fft_result) / N * 2  # 归一化并补偿窗口
    
    # 创建频率轴
    freq = (0:N-1) .* fs / N
    
    # 只保留正频率部分 (0到Nyquist频率)
    half_idx = Int(floor(N/2))
    freq = freq[1:half_idx]
    fft_magnitude = fft_magnitude[1:half_idx]
    
    return freq, fft_magnitude
end

# 5. 主函数：绘制信号和频谱
function plot_all_signals_and_spectrum()
    # 生成原始信号
    input_signal = generate_square_wave(t)
    
    # 生成曲顶抽样信号
    pulses = generate_natural_sampling_pulses(t, fs, input_signal, duty_cycle)
    
    # 计算频谱
    freq, spectrum = compute_spectrum(pulses, fs)
    
    # 创建图形（3个子图）
    fig = figure(figsize=(12, 10))
    
    # 子图1：原始信号（只显示前2个周期）
    subplot(3, 1, 1)
    display_range = findall(t .< 2 * T_signal)
    plot(t[display_range], input_signal[display_range], "b-", linewidth=1.5)
    title("原始信号 ($(signal_freq)Hz)")
    ylabel("幅度")
    grid(true)
    
    # 子图2：曲顶抽样信号（只显示前2个周期）
    subplot(3, 1, 2)
    plot(t[display_range], pulses[display_range], "r-", linewidth=1.0)
    title("曲顶抽样信号 (采样频率: $(fs)Hz, 占空比: $(round(duty_cycle*100, digits=1))%)")
    ylabel("幅度")
    grid(true)
    
    # 子图3：频谱
    subplot(3, 1, 3)
    plot(freq, spectrum, "g-", linewidth=1.5)
    # title("曲顶抽样信号频谱")
    # xlabel("频率 (Hz)")
    # ylabel("幅度")
    xlim(0, 5000)
    grid(true)
    
    # 添加总体标题
    

    return fig
end

# 运行绘图
fig = plot_all_signals_and_spectrum()


