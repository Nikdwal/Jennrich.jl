using Jennrich, Test, Random, Manifolds


@testset "CPD" begin
	Random.seed!(1)
	r = 3
	ð¯ = reshape(sum([kron(randn(9), randn(8), randn(7)) for i = 1:r]), (7,8,9))
	A, B, C = jennrich(ð¯, r)
	@test vec(ð¯) â sum([kron(C[:,i], B[:,i], A[:,i]) for i = 1:r])
end


@testset "(Láµ£,Láµ£,1)-BTD" begin
	Random.seed!(1)
	n = (9,8,7)
	L = [2,3]
	â³ = map(l -> Tucker(n, (l,l,1)), L)
	ð¯ = sum([embed(â³[i], TuckerPoint(randn(n...), (L[i],L[i],1))) for i = 1:length(L)])
	A,B,C = jennrich(ð¯, L)
	@test size(A) == (n[1], sum(L))
	@test size(B) == (n[2], sum(L))
	@test size(C) == (n[3], length(L))
	@test reshape(ð¯, n[1], :)' â [kron(C[:,1], B[:,1:L[1]]) kron(C[:,2], B[:,L[1]+1:end])] * A'
end
