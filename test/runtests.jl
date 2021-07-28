using Jennrich, Test, Random, Manifolds

@testset "(Lᵣ,Lᵣ,1)-BTD" begin
	Random.seed!(1)
	n = (9,8,7)
	L = [2,3]
	ℳ = map(l -> Tucker(n, (l,l,1)), L)
	𝔗 = sum([embed(ℳ[i], TuckerPoint(randn(n...), (L[i],L[i],1))) for i = 1:length(L)])
	A,B,C = jennrich(𝔗, L)
	@test size(A) == (n[1], sum(L))
	@test size(B) == (n[2], sum(L))
	@test size(C) == (n[3], length(L))
	@test reshape(𝔗, n[1], :)' ≈ [kron(C[:,1], B[:,1:L[1]]) kron(C[:,2], B[:,L[1]+1:end])] * A'
end
