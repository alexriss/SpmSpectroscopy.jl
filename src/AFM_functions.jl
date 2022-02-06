"""
    function deconvolve_sadar_jarvis(z::T, Δf::T, f₀::Float64, A::Float64, k::Float64)::Tuple{T,T} where T<:AbstractVector{Float64}
    
AFM force deconvolution using the Sadar-Jarvis method, as described in [Appl. Phys. Lett. 84, 1801 (2004)](https://aip.scitation.org/doi/10.1063/1.1667267).
Values for the tip-height `z` should be in ascending order. The corresponding values for the frequency shift are given in the vector `Δf`, along with
the experimental parameters `f₀` (resonance frequency), `A` (oscillation amplitude), and `k` (cantilever stiffness).

Based on MATLAB code from [Beilstein J. Nanotechnol. 3, 238 (2012)](https://www.beilstein-journals.org/bjnano/articles/3/27).
"""
function deconvolve_sader_jarvis(z::T, Δf::T, f₀::Float64, A::Float64, k::Float64)::Tuple{T,T} where T<:AbstractVector{Float64}
    @assert issorted(z)
    
    ω = Δf / f₀
    dω_dz = diff(ω) ./ diff(z)
    
    F = Vector{Float64}(undef, length(z) - 3)

    # we use `end-1` for both z and ω because dω_dz is shorter by one element
    for j=1:length(z)-3
        
        # adjust length of z, ω and dω_dz
        t = @view z[j+1:end-1]
        ω_ = @view ω[j+1:end-1]
        dω_dz_ = @view dω_dz[j+1:end]     
        
        integral = trapz(
            t,
            @. (1 + sqrt(A) / (8 * sqrt(π * (t - z[j])))) * ω_ - A^(3/2) / sqrt(2 * (t - z[j])) * dω_dz_
        )
                    
        # correction terms
        corr1 = ω[j] * (z[j+1] -z[j])
        corr2 = 2 * (sqrt(A) / (8 * sqrt(pi))) * ω[j] * sqrt(z[j+1] - z[j])
        corr3 = (-2) * (sqrt(A)^3 / sqrt(2)) * dω_dz[j] * sqrt(z[j+1] - z[j])
        F[j] = 2 * k * (corr1 + corr2 + corr3 + integral)
    end
    
    # adjust length of z
    z_ = z[1:length(F)]

    return z_, F
end


"""
    function deconvolve_matrix(z::T, Δf::T, f₀::Float64, A::Float64, k::Float64)::Tuple{T,T} where T<:AbstractVector{Float64}
    
AFM force deconvolution using the Matrix method, as described in [Appl. Phys. Lett. 78, 123 (2001)](https://aip.scitation.org/doi/10.1063/1.1335546).
Values for the tip-height `z` should be equally spaced and in ascending order. The corresponding values for the frequency shift are given in the vector `Δf`, along with
the experimental parameters `f₀` (resonance frequency), `A` (oscillation amplitude), and `k` (cantilever stiffness).
    
Based on MATLAB code from [Beilstein J. Nanotechnol. 3, 238 (2012)](https://www.beilstein-journals.org/bjnano/articles/3/27).
"""
function deconvolve_matrix(z::T, Δf::T, f₀::Float64, A::Float64, k::Float64)::Tuple{T,T} where T<:AbstractVector{Float64}
    @assert issorted(z)

    Δf_r = @view Δf[end:-1:1]
    Δz = z[2] - z[1]

    α = round(Int, A / Δz)

    N = length(z)
    # W = Array{Float64}(undef, N, N)
    W = zeros(Float64, N, N)

    for i = 1:N
        x = max(i - 2α, 1)
        for j = x:i
            W[i,j] = (f₀ / 2 / k) * (2 / π / A) * 2 / (2α + 1) * (sqrt((2α + 1) * (i - j + 1) - (i - j + 1)^2)-sqrt((2α + 1) * (i - j) - (i - j)^2))

            # τ1 = 1 - 2(i - j + 1) / (2α + 1)
            # τ2 = 1 - 2(i - j) / (2α + 1)
            # W[i,j] = (f₀ / 2 / k) * (2 / π / A) * (-sqrt(1 - τ2^2) + sqrt(1 - τ1^2))
        end
    end

    F = W \ Δf_r
    reverse!(F)
    return z[1:end], F  # return copy of z
end


"""
    function inflection_point_test(z::T, F::T, A::Float64; window_size::T2=nothing, polynomial_order::T2=nothing,
    window_size_deriv::T2=nothing, polynomial_order_deriv::T2=nothing)::NamedTuple where {T<:AbstractVector{Float64}, T2<:Union{Int,Nothing}}
    
Inflection point test on whether a AFM force deconvolution is reliable. The test has been described in [Nature Nanotechnology volume 13, 1088– (2018)](https://www.nature.com/articles/s41565-018-0277-x).
Values for the tip-height `z` should be equally spaced and in ascending order. The corresponding values for force (after deconvolution) are given in the vector `f`. The value for the oscillation amplitude
is given as `A`.

Window-size (`window_size`, should be an odd number) and polynomial order (`polynomial_order`) for the Savitzky-Golay smoothing will be automatically estimated if not given.
Similarly, the window-size (`window_size_deriv`, should be an odd number) and polynomial order (`polynomial_order_deriv`) for the Savitzky-Golay smoothing of the derivatives will be automatically estimated if not given.

Returns a named tuple:
  - `z_inp`: the z-values of the inflection points
  - `i_ifp`: the indices of the inflection points
  - `well_behaved`: characterization of each inflection point, i.e. 1 for well behaved passing the S_F test, 2 for passing the amplitude range test, and 0 for failing both tests.
  - `S_F`: the S_F value.
  - `amplitude_range`: amplitude ranges that possibly cause ill-posedness
  - `plot`: a plot with the fit, inflection points, as well as normalized and smoothed force and its derivatives (up to third order).

This function is highly experimental. It is highly dependent on the smoothing parameters.

Based on MATLAB code from [Journal of Applied Physics 127, 184301 (2020)](https://aip.scitation.org/doi/10.1063/5.0003291). As opposed to the spline fit suggested in the reference, Savitzky_Golay filtering is used here.

"""
function inflection_point_test(z::T, F::T, A::Float64; window_size::T2=nothing, polynomial_order::T2=nothing,
    window_size_deriv::T2=nothing, polynomial_order_deriv::T2=nothing)::NamedTuple where {T<:AbstractVector{Float64}, T2<:Union{Int,Nothing}}

    @assert length(z) == length(F)

    Δz = z[2] - z[1]

    if polynomial_order === nothing
        polynomial_order = 3
    end
    if polynomial_order_deriv === nothing
        polynomial_order_deriv = 3
    end
    if window_size === nothing
        window_size = round(Int, 80e-12 / Δz)
        if window_size % 2 == 0
            window_size += 1
        end
    end
    
    if window_size_deriv === nothing
        window_size_deriv = round(Int, 80e-12 / Δz)
        if window_size_deriv % 2 == 0
            window_size_deriv += 1
        end
    end
    
    
    # smoothed F, and its derivatives
    F_sg = savitzky_golay_filter(F, window_size, polynomial_order)
    F′_sg = savitzky_golay_filter(F, window_size, polynomial_order, deriv_order=1) ./ Δz
    F′_sg_sg = savitzky_golay_filter(F′_sg, window_size_deriv, polynomial_order_deriv)
    F′′_sg = savitzky_golay_filter(F, window_size, polynomial_order, deriv_order=2) ./ Δz^2
    F′′_sg_sg = savitzky_golay_filter(F′′_sg, window_size_deriv, polynomial_order_deriv)
    F′′′_sg = savitzky_golay_filter(F, window_size, polynomial_order, deriv_order=3)  ./ Δz^3
    F′′′_sg_sg = savitzky_golay_filter(F′′′_sg, window_size_deriv, polynomial_order_deriv)
    
    # limits from second last point to absolute max or min of F′_sg_sg
    N = length(F_sg)
    i_F′_min = findmin(F′_sg_sg)[2]
    i_F′_max = findmax(F′_sg_sg)[2]
    i_F′_extr = min(i_F′_min, i_F′_max)
    i_lim = [max(i_F′_extr, 2), N-1]  # cut outermost points
    
    # get inflection points by finding sign changes of second derivative
    i_ifp = Int[]
    for i in i_lim[2] : -1 : i_lim[1]  # from far to low
        if sign(F′′_sg_sg[i]) != sign(F′′_sg_sg[i+1])
            push!(i_ifp, i)
        end
    end
    
    z_ifp = z[i_ifp]
    z_ifp_half = z[i_ifp] ./ 2
    L_ifp = @. sqrt(abs.(-F′_sg_sg[i_ifp] / F′′′_sg_sg[i_ifp]))  # we have to put and `abs` here
    S_F = @. z[i_ifp]^2 / 4 * F′′′_sg_sg[i_ifp] ./ F′_sg_sg[i_ifp]
    
    
    well_behaved = map(enumerate(i_ifp)) do (i, _)
        if S_F[i] >= -1.
            return 1
        else
            z_str = @sprintf "%.3f nm"  z_ifp[i] * 1e9
            A_lower_str = @sprintf "%.1f pm"  L_ifp[i] * 1e12 
            A_upper_str = @sprintf "%.1f pm"  z_ifp_half[i] * 1e12 
            if A >= L_ifp[i] && A <= z_ifp_half[i]
                println("Inflection point at z = $(z_str): Amplitude inside of range possibly causing ill-posedness $(A_lower_str) to $(A_upper_str) ")
                return 0
            else
                println("Inflection point at z = $(z_str): Amplitude outside of range possibly causing ill-posedness $(A_lower_str) to $(A_upper_str) ")
                return 2
            end
        end
    end
    
    colors = map(well_behaved) do w
        if w == 1
            return "#20f060"
        elseif w == 2
            return "#90a020"
        elseif w == 0
            return "#f03060"
        else
            return "#a0a0a0"
        end
    end
    
    p1 = plot(z, F, label = "F", legend=:topleft)
    plot!(z, F_sg, label = "F_sg")
    scatter!(z[i_ifp], F_sg[i_ifp], label = "ifp", linewidth=0, ms=5, color=colors)
    
    p2 = plot(z, F_sg ./ maximum(abs.(F_sg)), label = "F_sg norm")
    plot!(z, F′_sg_sg ./ maximum(abs.(F′_sg_sg)), label = "F′_sg norm")
    plot!(z, F′′_sg_sg ./ maximum(abs.(F′′_sg_sg)), label = "F′′_sg norm")
    # plot!(z, F′′′_sg ./ maximum(abs.(F′′′_sg)), label = "F′′′_sg_sg norm")
    plot!(z, F′′′_sg_sg ./ maximum(abs.(F′′′_sg_sg)), label = "F′′′_sg_sg norm")
    
    p_all = plot(p1, p2, layout=(2,1))

    return (z_inp=z_ifp, i_ifp=i_ifp, well_behaved=well_behaved, S_F=S_F, amplitude_range=collect(zip(L_ifp, z_ifp_half)), plot=p_all)
end