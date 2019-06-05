using Strided
using BenchmarkTools
A = randn(4*256,4*256);
B = similar(A);
@btime $B .= ($A .* $A);
@btime @strided $B .= ($A .* $A);     
