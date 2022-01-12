using SpmSpectroscopy
using Test

@testset "loading" begin
    s = load_spectrum("Bias-Spectroscopy003.dat")
    @test all(s.position .≈ [ 1.1009e-8, -2.08172e-7, -9.98825e-9 ])
    @test s.header["Oscillation Control>Cut off frq (Hz)"] == "389"
    @test s.channel_names == [ "Bias calc", "Current", "Bias", "Phase", "Amplitude", "Frequency Shift", "Excitation"]
    @test s.channel_units == [ "V", "A", "V", "deg", "m", "Hz", "V"]
    @test s.z_feedback == false

    s = load_spectrum("Z-Spectroscopy__012.dat")
    @test all(s.position .≈ [  -1.99469e-8, -1.72054e-7, -1.88797e-8 ])
    @test s.channel_names == [ "Z rel", "Current", "Applied Voltage measured", "Bias", "X", "Y", "Z", "Phase", "Amplitude", "Frequency Shift", "Excitation", "Current [bwd]", "Applied Voltage measured [bwd]", "Bias [bwd]", "X [bwd]", "Y [bwd]", "Z [bwd]", "Phase [bwd]", "Amplitude [bwd]", "Frequency Shift [bwd]", "Excitation [bwd]" ]
    @test s.channel_units == [ "m", "A", "V", "V", "m", "m", "m", "deg", "m", "Hz", "V", "A", "V", "V", "m", "m", "m", "deg", "m", "Hz", "V" ]
end
