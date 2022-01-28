"""
    function trapz(X::AbstractVector{Float64}, Y::AbstractVector{Float64})::Float64

Trapezoidal integration over a function given by dicrete points in the arrays `Y` vs `X`,
where the spacing between the points in `X` is not necessarily constant.
"""
function trapz(X::AbstractVector{Float64}, Y::AbstractVector{Float64})::Float64
    @assert length(X) == length(Y)
    res = 0.
    for i = 1:length(X) - 1
        res += (X[i+1] - X[i]) * (Y[i] + Y[i+1])
    end
    return res / 2.
end


"""
    function rolling_mean(arr::AbstractArray{Float64}, n::Int64)::AbstractArray{Float64}

Computes the rolling mean over `n` points. The length of the output array is decreased by `n-1` points.

Adapted from:
https://stackoverflow.com/questions/59562325/moving-average-in-julia
"""
function rolling_mean(arr::AbstractArray{Float64}, n::Int64)::AbstractArray{Float64}
    so_far = sum(@view arr[1:n])
    out = zeros(Float64, length(arr) - n + 1)
    out[1] = so_far / n
    for (i, (start, stop)) in enumerate(zip(arr, @view arr[n+1:end]))
        so_far += stop - start
        out[i+1] = so_far / n
    end
    return out
end


"""
	function convolve_1d(u::AbstractVector, v::Vector)

One-dimensional convolution between the vectors `u` (size m) and `v` (size n). Returns a vector of size m+n-1.
"""
function convolve_1d(u::AbstractVector, v::AbstractVector)
    m = length(u)
    n = length(v)
    w = zeros(m + n - 1)
    @inbounds for j in 1:m, k in 1:n
        w[j+k-1] += u[j]*v[k]
    end
    #return w[n:end-n+1]
    return w
end


"""
	savitzky_golay_filter(y::AbstractVector, window_size::Integer, polynomial_order::Integer; deriv_order::Integer = 0, boundary_mode = :interpolation)

Apply Savitzky-Golay polynomial smoothing to input data `y` using a polynomial of order `polynomial_order` fit to a moving window `window_size` points wide. Optionally derivatives can be taken by specifying the `deriv_order` and the caller is responsible for the appropriate scaling by the point spacing. Handling of data within half the window size is specified by `boundary_mode`. When set to `:interpolation` the polynomial fit will be used; when set to `:nearest` the data will be padded using the edge values before convolution and the valid portion will then be returned.  
# References
1. Savitzky, A., & Golay, M. J. E. (1964). Smoothing and Differentiation of Data by Simplified Least Squares Procedures. Analytical Chemistry, 36(8), 1627–1639. https://doi.org/10.1021/ac60214a047
2. Steinier, J., Termonia, Y., & Deltour, J. (1972). Comments on Smoothing and differentiation of data by simplified least square procedure. Analytical Chemistry, 44(11), 1906–1909. https://doi.org/10.1021/ac60319a045
3. Press, W. H., & Teukolsky, S. A. (1990). Savitzky-Golay Smoothing Filters. Computers in Physics, 4(6), 669. https://doi.org/10.1063/1.4822961
"""
function savitzky_golay_filter(y::AbstractVector, window_size::Integer, polynomial_order::Integer; deriv_order::Integer = 0, boundary_mode = :interpolation)
    # from QLab.jl https://github.com/BBN-Q/Qlab.jl, Apache Software License

	# input validity checks
   	@assert isodd(window_size) "Window size must be an odd integer, i.e. fitting 2m + 1 points around the current value."
	@assert polynomial_order < window_size "Polynomial order must be less than the window size."
	@assert boundary_mode in (:interpolation, :nearest) "boundary_mode must be one of :interpolation, :nearest"

	# window size is 2m + 1 points
   	m = (window_size - 1) ÷ 2
	
	# build the Vandermonde design matrix A. Each row corresponds to a point in the fitting window -m:m
	# and each columns correspond to powers in the range 0:polynomial_order
   	fitting_points = -m:m
   	A = Matrix{Float64}(undef, window_size, polynomial_order + 1)
   	for i in 1:window_size, j in 1:polynomial_order + 1
        A[i,j] = fitting_points[i]^(j - 1)
    end

	if boundary_mode == :interpolation
		# for interpolation we'll want the full pseudo-inverse so we can calculate all the fit values at the edges
		# Ap = y
		C = pinv(A)

		# the filter coefficients are the rows of `C`
		filter_coeffs = C[deriv_order + 1,:] * factorial(deriv_order)

		# convolve with the filter coefficients with a couple extra steps:
		# 1. because of convolution will reverse coefficients we flip before
		# 2. c = conv(a,b) will return a vector of length(c) = length(a) + length(b) - 1 so we chop off the first and last m points
		smoothed = convolve_1d(reverse(filter_coeffs), y)[m+1:end-m]

		# for interpolation edge handling calculate the full fits
		if deriv_order == 0
			# if we are just smoothing then we can use the design and coefficient matrix as is
			AC = A*C
			smoothed[1:m] = (AC*y[1:window_size])[1:m]
			smoothed[end-m+1:end] = (AC*y[end-window_size+1:end])[end-m+1:end]
		else
			# otherwise we need to differentiate the polynomial coefficients
			# first m points
			p = C * y[1:window_size]
			for _ in 1:deriv_order
				p = [(i-1)*p[i] for i in 2:length(p)]
			end
			smoothed[1:m] = A[1:m, 1:size(A,2)-deriv_order]*p
			# last m points
			p = C * y[end-window_size+1:end]
			for _ in 1:deriv_order
				p = [(i-1)*p[i] for i in 2:length(p)]
			end
			smoothed[end-m+1:end] = A[m+2:end, 1:size(A,2)-deriv_order]*p

		end

		return smoothed

	elseif boundary_mode == :nearest

		# here we only need a single set of coefficients and so we can least-squares solve AᵀCᵀ = I for for a single row by picking a single column out of I
		Icol = zeros(Float64, polynomial_order+1, 1)
		Icol[deriv_order + 1] = 1.0
		filter_coeffs = transpose(A) \ Icol

		# pad the signal with the endpoints
		padded_y = [y[1] * ones(m); vec(y); y[end] * ones(m)]

		# convolve with filter
		smoothed = convolve_1d(filter_coeffs[end:-1:1], padded_y)

		# and return the valid midsection
		return smoothed[window_size:end-2*m]

	end

end