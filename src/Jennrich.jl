module Jennrich
using Manifolds, LinearAlgebra

export jennrich

"""
	jennrich(A, r :: Int)

Compute the rank r canonical polyadic decomposition of A using a pencil-based algorithm
"""
jennrich(A, r :: Int) = jennrich(A, ones(Int, r))

"""
	jennrich(A, L :: AbstractVector{Int})

Compute the (Lᵣ,Lᵣ,1)-BTD of A using a pencil-based algorithm
"""
function jennrich(𝔗 :: AbstractArray{<:Real, 3}, L :: AbstractVector{Int})
	n₁, n₂, n₃ = size(𝔗)
	@assert n₁ ≥ sum(L) && n₂ ≥ sum(L)

	𝔗Tuck = TuckerPoint(𝔗, (sum(L), sum(L), 2))
	𝔊     = 𝔗Tuck.hosvd.core

	evd   = eigen(𝔊[:,:,1] * 𝔊[:,:,2]^-1)
	makeReal!(evd)
	A :: Matrix{eltype(𝔗)} = 𝔗Tuck.hosvd.U[1] * real(evd.vectors)

	B⨀C   = (A \ reshape(𝔗, n₁, :))'
	B     = Matrix{eltype(A)}(undef, n₂, sum(L))
	C     = Matrix{eltype(A)}(undef, n₃, length(L))
	blkInd = cumsum(vcat(1, L))
	for i ∈ 1:length(L)
		# Bᵢ⨀Cᵢ is the concatenation of rank 1 matrices Bᵢ[:,j] ⊗ cᵢ
		# Compute the rank 1 SVD of these matrices (they all have cᵢ as a right singular vector)
		from, to = blkInd[i], blkInd[i+1] - 1
		Bᵢ⨀Cᵢ = B⨀C[:,from:to] 
		rk1mtx = reshape(Bᵢ⨀Cᵢ * randn(L[i]), (n₂, n₃))
		C[:,i] = leftSingularVector(rk1mtx', maxiter=2) # right singular vector
		for j = 1:L[i]
			rk1mtx = reshape(Bᵢ⨀Cᵢ[:,j], (n₂, n₃))
			B[:,from-1+j] = rk1mtx * C[:,i]  # left singular vector * singular value
		end
	end
	A, B, C
end

# Compute the dominant left singular vector of A
# This only works well if the second singular value of A is significantly less than the first
function leftSingularVector(A; maxiter=1)
	σv = randn(size(A,2))
	u = normalize(A * σv)
	for i = 1:maxiter - 1
		mul!(σv, A', u)
		normalize!(mul!(u, A, σv))
	end
	u
end

function makeReal!(evd::Union{GeneralizedEigen, Eigen})
	λ, V = evd
	i = firstindex(λ)
	while i ≤ lastindex(λ)
		if imag(λ[i]) ≠ 0
			λ[i+1] = conj(λ[i+1])
			V[:,i] = real.(V[:,i])
			V[:,i+1] = imag.(V[:,i+1])
			i += 1
		end
		i += 1
	end
	evd
end

end # module
