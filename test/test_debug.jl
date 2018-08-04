# using BenchmarkTools

# @inline bis(a,b) = (a+b)/2
# # Note, the benchmark that follows uses atol = sqrt(eps())
# # but I decided to update it later to reflect the value of the upper bracket.
# function bisection(f, a_, b_, atol = 2eps(promote_type(typeof(b_),Float64)(b_)))
#     c = bis(a_,b_)
#     z = f(c)
#     if z > 0 #
#         b = c
#         a = typeof(b)(a_)
#     else
#         a = c
#         b = typeof(a)(b_)
#     end
#     while abs(a - b) > atol
#         c = bis(a,b)
#         if f(c) > 0 #
#             b = c
#         else
#             a = c
#         end
#     end
#     a, b
# end

# f(x) = exp(x) - x^4

# bisection(f,8,9)

# function bisect1(fun, xl::Real, xu::Real, tolx::Real, tolf::Real)
#     if (xl > xu)
#         xl, xu = xu, xl
#     end
#     fl = fun(xl);
#     fu = fun(xu);
#     @assert fl*fu < 0.0 "Need to get a bracket"
#     if fl == 0.0
#         return xl, xl;
#     end
#     if fu == 0.0
#         return xu, xu;
#     end
#     while true
#         xr = (xu + xl) / 2 ; # bisect interval
#         fr = fun(xr); # value at the midpoint
#         if fr * fl < 0 # (fr < 0.0 && fl > 0.0) || (fr > 0.0 && fl < 0.0)
#             xu, fu = xr, fr;# upper --> midpoint
#         else
#             xl, fl = xr, fr;# lower --> midpoint
#         end
#         (abs(fr) <= tolf) && return xl, xu;
#         ((xu-xl) < tolx) && return xl, xu;
#     end
#     return xl, xu;
# end

# tolx = 2eps(9.0)
# tolf = √eps()
# bisect1(f, 8.0, 9.0, tolx, tolf)

# """
#     bisection(f, a, b; fa = f(a), fb = f(b), ftol, wtol)

# Bisection algorithm for finding the root ``f(x) ≈ 0`` within the initial bracket
# `[a,b]`.

# Returns a named tuple

# `(x = x, fx = f(x), isroot = ::Bool, iter = ::Int, ismaxiter = ::Bool)`.

# Terminates when either

# 1. `abs(f(x)) < ftol` (`isroot = true`),
# 2. the width of the bracket is `≤wtol` (`isroot = false`),
# 3. `maxiter` number of iterations is reached. (`isroot = false, maxiter = true`).

# which are tested for in the above order. Therefore, care should be taken not to make `wtol` too large.

# """
# function bisect2(f, a::Real, b::Real; fa::Real = f(a), fb::Real = f(b),
#                    ftol = √eps(), wtol = 0, maxiter = 100)
#     @assert fa * fb ≤ 0 "initial values don't bracket zero"
#     @assert isfinite(a) && isfinite(b)
#     _bisection(f, float.(promote(a, b, fa, fb, ftol, wtol))..., maxiter)
# end

# function _bisection(f, a, b, fa, fb, ftol, wtol, maxiter)
#     iter = 0
#     abs(fa) < ftol && return (x = a, fx = fa, isroot = true, iter = iter, ismaxiter = false)
#     abs(fb) < ftol && return (x = b, fx = fb, isroot = true, iter = iter, ismaxiter = false)
#     while true
#         iter += 1
#         m = (a + b) / 2
#         fm = f(m)
#         abs(fm) < ftol && return (x = m, fx = fm, isroot = true, iter = iter, ismaxiter = false)
#         abs(b-a) ≤ wtol && return (x = m, fx = fm, isroot = false, iter = iter, ismaxiter = false)
#         if fa * fm > 0
#             a, fa = m, fm
#         else
#             b, fb = m, fm
#         end
#         iter == maxiter && return (x = m, fx = fm, isroot = false, iter = iter, ismaxiter = true)
#     end
# end


