
module mtfieldtm1
using FinEtools
using Test
using BenchmarkTools
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 20, 30, 40;

    fens,fes  = H8block(L,W,t, nl,nw,nt)
    @show count(fens)
    geom  =  NodalField(fens.xyz)

@btime (l1 = selectnode($fens; box=[0.0 $L 0.0 $W 0.0 $t], inflate = $t/100.0))


end
end
using Main.mtfieldtm1
mtfieldtm1.test()
#
# module m1
# function test()
#     @inline inrange(rangelo,rangehi,x) = ((x>=rangelo) && (x<=rangehi));
#     N = 26_000
#     sdim = 3
#     v = rand(N, 3)
#     abox = vec(Float64[0.0 0.0 0.0 1.0 0.0 1.0])
#     nn = 1
#     @time begin
#         for j in 1:size(v,1)
#             match = true
#             for i in 1:sdim
#                 if (!inrange(abox[2*i-1],abox[2*i], v[j, i]))
#                     match =  false; break
#                 end
#             end
#             if match
#                 nn = nn + 1;
#             end
#         end
#     end
# end
# end
# using .m1
# m1.test()

# module m2
#  using FinEtools
#  using Test
#  function vselect(v::Matrix{Float64}; kwargs...)
#
#      # Helper functions
#      @inline inrange(rangelo,rangehi,x) = ((x>=rangelo) && (x<=rangehi));
#
#      # Extract arguments
#      box = nothing; inflate = 0.0;
#      for apair in pairs(kwargs)
#          sy, val = apair
#          if sy == :box
#              box = val
#          elseif sy == :inflate
#              inflate = val
#          end
#      end
#
#      # Did we get an inflate value
#      inflatevalue =0.0;
#      if inflate != nothing
#          inflatevalue = inflate;
#      end
#
#      # Initialize the output list
#      vlist= zeros(FInt, size(v,1)); nn= 0;
#
#
#      # Process the different options
#      if box != nothing
#          sdim = size(v,2)
#          dim = Int(round(length(box)/2.))::FInt;
#          @assert dim == sdim "Dimension of box not matched to dimension of array of vertices"
#          abox = FFltVec(vec(box))
#          inflatebox!(abox, inflatevalue)
#          for j in 1:size(v,1)
#              match = true
#              for i in 1:sdim
#                  loc = v[j, i]
#                  if (!inrange(abox[2*i-1],abox[2*i],loc))
#                      match =  false; break
#                  end
#              end
#              if match
#                  nn = nn + 1; vlist[nn] = j;
#              end
#          end
#      end
#      if (nn==0)
#          vlist = FInt[];# nothing matched
#      else
#          vlist = vlist[1:nn];
#      end
#      return vlist
#  end
#  function test()
#      N = 26_000
#      sdim = 3
#      v = rand(N, 3)
#      abox = vec(Float64[0.0 0.0 0.0 1.0 0.0 1.0])
#      @time vselect(v; box = abox)
#  end
#  end
#  using .m2
#  m2.test()
