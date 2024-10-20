################## RF Analysis #######################

using FFTW

"""
    complex_signal(cf::ControlField)

Calculates complex RF signal in Hz

"""
complex_signal(cf::ControlField) = vec(cf.B1x + im * cf.B1y)

"""
    max_peak_amp(B1::Vector{ComplexF64})

Returns the maximum RF amplitude in Hz
"""
max_peak_amp(B1::Vector{ComplexF64}) = round(maximum(abs.(B1)), digits = 4)

"""
    max_peak_amp_tesla(B1::Vector{ComplexF64})

Returns the maximum RF amplitude in Tesla
"""
max_peak_amp_tesla(B1::Vector{ComplexF64}) = round((maximum(abs.(B1))/γ_¹H)*1e6, digits = 4)


"""
    amplitudes_and_phases(cf::ControlField)

Calculates amplitudes and phases in Hz and radians

"""
function amplitudes_and_phases(cf::ControlField)
    B1 = complex_signal(cf)
    return abs.(B1), angle.(B1)
end


#########################################################
################## Pulse Analysis #######################
#########################################################
"""
    RF_pulse_analysis(cf::ControlField; attenuation_ref = 0.0, B1_ref = 1.0, power_ref = 500.0)

Calculate calibration analysis of shaped pulse

    ### Parameters:
    - `cf::ControlField`: The control field object containing the RF waveform.
    - `attenuation_ref::Float64`: Reference attenuation in dB (default is 0.0 dB).
    - `B1_ref::Float64`: Reference RF field strength in Tesla (default is 1.0 T).
    - `power_ref::Float64`: Reference RF power in Watts (default is 500.0 W).
    
    ### Returns:
    A tuple of calculated values:
    - `max_amp`: Maximum RF amplitude in Hertz (Hz).
    - `max_amp_tesla`: Maximum RF amplitude in microtesla (µT).
    - `attenuation_B1`: Attenuation in decibels (dB).
    - `power_max_B1`: Maximum power in Watts (W).
    - `power_average`: Average power in Watts (W).
    - `pulse_energy`: Pulse energy in Joules (J).

"""
function RF_pulse_analysis(cf::ControlField; attenuation_ref = 0.0, B1_ref = 1.0, power_ref = 500.0)
    B1 = complex_signal(cf)
    println("---------- RF Analysis ---------- \n")
    # RF maximum amplitude
    max_amp = max_peak_amp(B1)
    println("Pulse Peak Amplitude = $max_amp [Hz]")
    max_amp_tesla = max_peak_amp_tesla(B1)
    println("Pulse Peak Amplitude = $max_amp_tesla [μT]")

    # Attenuation
    ratio_B1 = B1_ref/max_amp_tesla # [T]
    attenuation_B1 = round(attenuation_ref + 20*log(ratio_B1), digits = 4)
    println("Attenuation corresponding to maximum amplitude: $attenuation_B1 [dB]")

    # Power
    power_max_B1 = round(power_ref/(10^(attenuation_B1/20)), digits = 4)
    println("Maximum power = $power_max_B1 [W]")
    instantaneous_power_B1 = (abs.(B1).^2) / maximum(abs.(B1).^2) * power_max_B1
    power_average = round(1/length(B1)*sum(instantaneous_power_B1), digits = 4)
    println("Average power = $power_average [W]")

    # Energy
    pulse_energy = round(cf.t_control*power_average, digits = 4) # W = J/s
    println("Pulse energy = $pulse_energy [J]")

    return max_amp, max_amp_tesla, attenuation_B1, power_max_B1, power_average, pulse_energy
end


#########################################################
################# Fourier Analysis ######################
#########################################################
# function fast_fourier_transform(cf::ControlField)
#     # RF Pulse 
#     N = length(cf.B1x)
#     B1 = complex_signal(cf)
#     time = range(0.0, stop=cf.t_control, length=N)

#     # Compute FFT
#     Δt = cf.t_control/N
#     sampling_rate = 1/Δt
#     fft_RF = FFTW.fft(B1)
#     frequencies = fftfreq(N, sampling_rate) |> fftshift

#     # Plot
#     plot_layout = @layout [a; b]
#     pTime = plot()
#         plot!(pTime, time, abs.(B1), label = false, ylabel = "Magnitude", xlabel = "Time [s]", title = "Time-Domain Pulse")
#     pFreq = plot()
#         plot!(pFreq, frequencies, fftshift(abs.(fft_RF)), label = false, ylabel = "Magnitude", xlabel = "Frequency [kHz]", title = "Frequency-Domain Pulse")
#     p = plot(pTime, pFreq, layout = plot_layout)

#     # cf = res_grape.control_field 
#     # pFFT = fast_fourier_transform(cf)

#     return p
# end


#########################################################
############## Bruker Data Conversion ###################
#########################################################
"""
    bruker_normalized_amplitudes_and_phases(cf::ControlField)

Calculates amplitudes and phases normalized to 100 and 180 deg for Bruker implementation on TopSpin.
All negative phase values are added a 360deg phase
    
    ### Parameters:
    - `cf::ControlField`: The control field object containing the RF waveform.
    
    ### Returns:
    A tuple of normalized amplitudes and phases:
    - `norm_amp`: Maximum RF amplitude in Hertz (Hz).
    - `norm_phase`: Normalize phases in degrees. 

"""
function bruker_normalized_amplitudes_and_phases(cf::ControlField)
    (ampli, phase) = amplitudes_and_phases(cf)
    norm_ampli = (ampli ./ maximum(ampli)) * 100
    norm_phase = phase .* (180/π) # Convert phase from radians to degree
    norm_phase = norm_phase .+ 360 .* (norm_phase .< 0) # Adding a 360° phase to negative phase values
    return norm_ampli, norm_phase
end


"""
    bandwidth_factor(cf::ControlField)

Excitation frequency range

"""
function bandwidth_factor()
    # RF Pulse 
    N  = length(cf.B1x)
    B1 = complex_signal(cf)
    Δt = cf.t_control/N
    t  = range(0.0, cf.t_control, length=N)

    # Perform FFT
end 

"""
    integral_factor(cf::ControlField)

Power or energy of the RF pulse integrated over its duration

"""
function integral_factor(cf::ControlField)
    B1 = abs.(complex_signal(cf))
    norm_amp = B1./maximum(B1)
    rf_time = range(0, cf.t_control, length(1:length(cf.B1x)))
    int_fac = integrate(rf_time, norm_amp)
    return round(int_fac, digits=4)
end

