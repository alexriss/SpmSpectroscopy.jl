using SpmSpectroscopy
using Test

@testset "loading" begin
    s = load_spectrum("Bias-Spectroscopy003.dat")
    @test all(s.position .≈ [ 1.1009e-8, -2.08172e-7, -9.98825e-9 ])
    @test s.header["Oscillation Control>Cut off frq (Hz)"] == "389"
    @test s.channel_names == [ "Bias calc", "Current", "Bias", "Phase", "Amplitude", "Frequency Shift", "Excitation"]
    @test s.channel_units == [ "V", "A", "V", "deg", "m", "Hz", "V"]
    @test s.z_feedback == false

    s = load_spectrum("Z-Spectroscopy__012.dat", index_column=true)
    @test all(s.position .≈ [  -1.99469e-8, -1.72054e-7, -1.88797e-8 ])
    @test s.channel_names == [ "Z rel", "Current", "Applied Voltage measured", "Bias", "X", "Y", "Z", "Phase", "Amplitude", "Frequency Shift", "Excitation", "Current [bwd]", "Applied Voltage measured [bwd]", "Bias [bwd]", "X [bwd]", "Y [bwd]", "Z [bwd]", "Phase [bwd]", "Amplitude [bwd]", "Frequency Shift [bwd]", "Excitation [bwd]", "Index" ]
    @test s.channel_units == [ "m", "A", "V", "V", "m", "m", "m", "deg", "m", "Hz", "V", "A", "V", "V", "m", "m", "m", "deg", "m", "Hz", "V", "" ]
    @test s.data[!,"Index"] == collect(1:size(s.data,1))

    s = load_spectrum("Z-Spectroscopy__012.dat", index_column=true, index_column_type=Float64)
    @test s.data[3,"Index"] === 3.0

    s = load_spectrum("Z-Spectroscopy002.dat", remove_missing=true)
    @test size(s.data) == (77,17)

    s = load_spectrum("Z-Spectroscopy002.dat")  # default for remove_missing is false
    @test size(s.data) == (256,17)

    # Base.show
    s = load_spectrum("Z-Spectroscopy__012.dat", index_column=true, index_column_type=Float64)
    io = IOBuffer()
    print(IOContext(io, :compact => false), s)
    @test String(take!(io)) == """SpmSpectrum("Z-Spectroscopy__012.dat", Experiment: "Z spectroscopy", 22 channels, 128 points)"""
    print(IOContext(io, :compact => true), s)
    @test String(take!(io)) == """SpmSpectrum("Z-Spectroscopy__012.dat")"""
    s = load_spectrum("Z-Spectroscopy__012.dat", index_column=true, index_column_type=Float64, header_only=true)
    io = IOBuffer()
    print(IOContext(io, :compact => false), s)
    @test String(take!(io)) == """SpmSpectrum("Z-Spectroscopy__012.dat", Experiment: "Z spectroscopy", 22 channels, 0 points)"""
    print(IOContext(io, :compact => true), s)
    @test String(take!(io)) == """SpmSpectrum("Z-Spectroscopy__012.dat")"""
end

@testset "background corrections" begin
    x = [1.,2.,3.]
    y = [4.,5.,6.]
    correct_background!(x, y, no_correction)
    @test x == [1.,2.,3.]
    @test y == [4.,5.,6.]
    correct_background!(x, y, subtract_minimum)
    @test x == [1.,2.,3.]
    @test all(y .≈ [0.,1.,2.])
    x = [1.,2.,3.]
    y = [4.,5.,7.]
    correct_background!(x, y, linear_fit)
    @test all(y .≈ [0.5,0.,0.5])

    # some special cases
    x = [0., 0.]
    y = [2., 2.]
    correct_background!(x, y, linear_fit)
    @test y ≈ [0., 0.]
    x = [0.]
    y = [0.]
    correct_background!(x, y, linear_fit)
    @test y ≈ [0.]
    x = Float64[]
    y = Float64[]
    correct_background!(x, y, linear_fit)
    @test length(y) == 0
end

@testset "spectrum manipulations" begin
    @test all(SpmSpectroscopy.rolling_mean([1.,2.,5.,6.,7.], 3) .≈ [1.0+2.0+5.0, 2.0+5.0+6.0, 5.0+6.0+7.0] ./3)
    @test SpmSpectroscopy.trapz([1.,2,3], [1.,1,1]) ≈ 2.0
    @test SpmSpectroscopy.trapz([1.,2,4], [1.,1,1]) ≈ 3.0

    @test SpmSpectroscopy.trapz([1.,2,3], [1.,1,1]) ≈ 2.0
    @test SpmSpectroscopy.trapz([1.,2,4], [1.,1,1]) ≈ 3.0

    @test all(SpmSpectroscopy.convolve_1d([1,4,6], [6,4,5]) .≈ [6.0, 28.0, 57.0, 44.0, 30.0])

    a = [74, 36, 47, 14, 3, 89, 84, 73, 46, 47, 43, 41, 40, 74, 63, 51, 38, 6, 19, 44, 7, 29, 23, 29, 90, 61, 40, 74, 67, 28, 67, 72, 54, 97, 7, 9, 30, 53, 54, 7, 62, 59, 45, 11, 24, 43, 5, 100, 38, 66]
    a_sg = [74.30952381, 43.16666667, 25.5952381 , 21.76190476, 38.76190476,
    58.0952381 , 75.47619048, 76.52380952, 55.38095238, 45.57142857,
    37.76190476, 44.28571429, 52.        , 59.52380952, 63.38095238,
    51.66666667, 29.42857143, 25.19047619, 19.42857143, 21.28571429,
    25.85714286, 18.28571429, 28.0952381 , 50.33333333, 54.9047619 ,
    63.61904762, 68.9047619 , 53.        , 54.0952381 , 59.52380952,
    51.9047619 , 69.38095238, 73.33333333, 52.0952381 , 32.71428571,
    24.71428571, 26.52380952, 37.38095238, 41.80952381, 44.33333333,
    47.57142857, 45.38095238, 42.52380952, 31.57142857, 15.42857143,
    30.57142857, 44.04761905, 54.88095238, 61.52380952, 61.30952381]  #using scipy.signal.savgol_filter(a, 7, 3, mode="interp")
    @test all(abs.(SpmSpectroscopy.savitzky_golay_filter(a, 7, 3) .- a_sg) .< 1e-6)

    a_sg = [ 13.4047619 ,  13.57142857,  13.73809524,  13.9047619,
         5.33333333,  -3.61904762, -12.30952381, -10.76190476,
         2.52380952,   3.92857143,   7.11904762,   3.14285714,
        -0.35714286,  -4.76190476,  -9.33333333,  -5.04761905,
         6.35714286,   3.69047619,   4.14285714,   1.21428571,
        -1.71428571,   8.07142857,   6.16666667,  -5.23809524,
        -2.73809524,  -4.38095238,  -6.66666667,   4.        ,
         2.16666667,  -1.04761905,   6.83333333,  -6.69047619,
       -12.80952381,  -2.04761905,   6.64285714,   9.35714286,
         5.0952381 ,  -2.83333333,  -1.33333333,  -0.02380952,
        -3.        ,  -3.97619048,  -3.33333333,   2.        ,
        12.78571429,   3.71428571,  -1.52380952,  -4.19047619,
        -6.85714286,  -9.52380952] #using scipy.signal.savgol_filter(a, 7, 3, deriv=2, mode="interp")
    @test all(abs.(SpmSpectroscopy.savitzky_golay_filter(a, 7, 3, deriv_order=2) .- a_sg) .< 1e-6)
end

@testset "AFM data manipulations" begin
    s = load_spectrum("Z-Spectroscopy__012.dat", index_column=true, index_column_type=Float64)
    z = s.data[:, "Z rel"]
    df = s.data[:, "Frequency Shift"]

    # deconvolutions

    # F_ml is the result using MATLAB code
    F_ml = 1e-9 .* [-0.1116,  -0.1125,  -0.1141,  -0.1190,  -0.1078,  -0.1168,  -0.1110,  -0.1248,  -0.1175,  -0.1166,  -0.1263,  -0.1278,  -0.1232,  -0.1280,  -0.1336,  -0.1310,  -0.1353,  -0.1314,  -0.1418,  -0.1308,  -0.1396,  -0.1375,  -0.1329,  -0.1427,  -0.1338,  -0.1355,  -0.1348,  -0.1383,  -0.1305,  -0.1327,  -0.1279,  -0.1240,  -0.1282,  -0.1234,  -0.1239,  -0.1166,  -0.1169,  -0.1158,  -0.1135,  -0.1103,  -0.1072,  -0.1122,  -0.1058,  -0.1033,  -0.1023,  -0.0981,  -0.0968,  -0.0958,  -0.0926,  -0.0910,  -0.0924,  -0.0831,  -0.0841,  -0.0832,  -0.0845,  -0.0733,  -0.0771,  -0.0810,  -0.0696,  -0.0772,  -0.0741,  -0.0666,  -0.0674,  -0.0664,  -0.0676,  -0.0642,  -0.0622,  -0.0644,  -0.0549,  -0.0613,  -0.0557,  -0.0572,  -0.0564,  -0.0444,  -0.0578,  -0.0505,  -0.0492,  -0.0537,  -0.0397,  -0.0510,  -0.0437,  -0.0425,  -0.0414,  -0.0450,  -0.0404,  -0.0352,  -0.0424,  -0.0377,  -0.0346,  -0.0320,  -0.0397,  -0.0299,  -0.0363,  -0.0294,  -0.0275,  -0.0288,  -0.0308,  -0.0272,  -0.0249,  -0.0239,  -0.0269,  -0.0225,  -0.0226,  -0.0164,  -0.0262,  -0.0145,  -0.0206,  -0.0195,  -0.0165,  -0.0185,  -0.0133,  -0.0145,  -0.0148,  -0.0162,  -0.0078,  -0.0151,  -0.0066,  -0.0129,  -0.0064,  -0.0053,  -0.0120,  -0.0011,  -0.0066,  -0.0066,   0.0017]
    z_, F = SpmSpectroscopy.deconvolve_sader_jarvis(z, df, 30000., 60e-12, 1800.)
    @test all(abs.(F .- F_ml) .< 1e-13)

    # padding
    z_, F = SpmSpectroscopy.deconvolve_sader_jarvis(z, df, 30000., 60e-12, 1800., pad=true, val=missing)
    @test all(abs.(skipmissing(F) .- F_ml) .< 1e-13)
    @test length(F) == length(z)
    @test all(ismissing.(F[end-2:end]))
    z_, F = SpmSpectroscopy.deconvolve_sader_jarvis(z, df, 30000., 60e-12, 1800., pad=true, val=0.7)
    @test length(F) == length(z)
    @test all(F[end-2:end] .≈ 0.7)
    z_, F = SpmSpectroscopy.deconvolve_sader_jarvis(z, df, 30000., 60e-12, 1800., pad=true, val=NaN)
    @test length(F) == length(z)
    @test all(isnan.(F[end-2:end]))

    # Matrix deconvolution
    F_ml = 1e-9 .* [-0.1303,  -0.1390,  -0.1271,  -0.1428,  -0.1313,  -0.1306,  -0.1337,  -0.1387,  -0.1313,  -0.1314,  -0.1344,  -0.1205,  -0.0841,  -0.1297,  -0.1427,  -0.1420,  -0.1539,  -0.1459,  -0.1520,  -0.1492,  -0.1527,  -0.1568,  -0.1546,  -0.1539,  -0.1500,  -0.1512,  -0.1540,  -0.1487,  -0.1495,  -0.1452,  -0.1363,  -0.1389,  -0.1380,  -0.1388,  -0.1334,  -0.1185,  -0.1274,  -0.1260,  -0.1188,  -0.1229,  -0.1092,  -0.1226,  -0.1103,  -0.1038,  -0.1067,  -0.0999,  -0.0970,  -0.0930,  -0.0867,  -0.0772,  -0.0680,  -0.1200,  -0.1145,  -0.1069,  -0.1078,  -0.0948,  -0.0951,  -0.0954,  -0.0902,  -0.0929,  -0.0858,  -0.0789,  -0.0815,  -0.0786,  -0.0776,  -0.0704,  -0.0747,  -0.0668,  -0.0641,  -0.0678,  -0.0628,  -0.0648,  -0.0577,  -0.0522,  -0.0623,  -0.0580,  -0.0531,  -0.0576,  -0.0454,  -0.0565,  -0.0465,  -0.0452,  -0.0503,  -0.0458,  -0.0453,  -0.0446,  -0.0445,  -0.0444,  -0.0532,  -0.1005,  -0.0660,  -0.0523,  -0.0512,  -0.0423,  -0.0393,  -0.0395,  -0.0389,  -0.0345,  -0.0319,  -0.0315,  -0.0323,  -0.0284,  -0.0271,  -0.0245,  -0.0292,  -0.0217,  -0.0264,  -0.0251,  -0.0235,  -0.0240,  -0.0207,  -0.0223,  -0.0232,  -0.0229,  -0.0186,  -0.0231,  -0.0187,  -0.0226,  -0.0184,  -0.0202,  -0.0240,  -0.0185,  -0.0245,  -0.0246,  -0.0224,  -0.0306,  -0.0402,  -0.0590]
    z_, F = SpmSpectroscopy.deconvolve_matrix(z, df, 30000., 60e-12, 1800.)
    @test all(abs.(F .- F_ml) .< 1e-13)
    @test length(F) == length(z)  # no padding required

    # we skip the inflection point test for now, it is so sensitive to smoothing
end
