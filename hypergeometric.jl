using Nemo
import Nemo: libarb, libflint, arf_struct, arb_struct
import Nemo: zero!, one!, add!, sub!, mul!, div!, coeff,
       setcoeff!, mullow, gamma, rgamma, solve, overlaps

const ARF_PREC_EXACT = typemax(Int)

# TODO the naked arb is the same as Nemo.acb_struct, maybe use that instead
mutable struct narb
  mid_exp::Int # fmpz
  mid_size::UInt # mp_size_t
  mid_d1::UInt # mantissa_struct
  mid_d2::UInt
  rad_exp::Int # fmpz
  rad_man::UInt

  function narb()
    z = new()
    ccall((:arb_init, libarb), Nothing,
        (Ref{narb}, ),
        z)
    finalizer(_arb_clear_fn, z)
    return z
  end
end

function _arb_clear_fn(x::narb)
  ccall((:arb_clear, libarb), Nothing, (Ref{narb}, ), x)
end

mutable struct nacb
  real_mid_exp::Int # fmpz
  real_mid_size::UInt # mp_size_t
  real_mid_d1::Int # mantissa_struct
  real_mid_d2::Int
  real_rad_exp::Int # fmpz
  real_rad_man::UInt
  imag_mid_exp::Int # fmpz
  imag_mid_size::UInt # mp_size_t
  imag_mid_d1::Int # mantissa_struct
  imag_mid_d2::Int
  imag_rad_exp::Int # fmpz
  imag_rad_man::UInt

  function nacb()
    z = new()
    ccall((:acb_init, libarb), Nothing,
        (Ref{nacb}, ),
        z)
    finalizer(_acb_clear_fn, z)
    return z
  end
end

function _acb_clear_fn(x::nacb)
  ccall((:acb_clear, libarb), Nothing, (Ref{nacb}, ), x)
end

mutable struct nacb_poly
  coeffs::Ptr{Nothing}
  length::Int
  alloc::Int

  function nacb_poly()
    z = new()
    ccall((:acb_poly_init, libarb), Nothing, (Ref{nacb_poly}, ), z)
    finalizer(_acb_poly_clear_fn, z)
    return z
  end
end

function _acb_poly_clear_fn(x::nacb_poly)
  ccall((:acb_poly_clear, libarb), Nothing, (Ref{nacb_poly}, ), x)
end


mutable struct nacb_mat
  entries::Ptr{Nothing}
  r::Int
  c::Int
  rows::Ptr{Nothing}

  function nacb_mat(r::Int, c::Int)
    z = new()
    ccall((:acb_mat_init, libarb), Nothing,
          (Ref{nacb_mat}, Int, Int),
          z, r, c)
    finalizer(_acb_mat_clear_fn, z)
    return z
  end
end

function _acb_mat_clear_fn(a::nacb_mat)
  ccall((:acb_mat_clear, libarb), Nothing, (Ref{nacb_mat}, ), a)
end

################################################################################

function nacb(x::Int)
  z = nacb()
  ccall((:acb_set_si, libarb), Nothing,
        (Ref{nacb}, Int),
        z, x)
  return z
end

function nacb(x::narb)
  z = nacb()
  ccall((:acb_set_arb, libarb), Nothing,
        (Ref{nacb}, Ref{narb}),
        z, x)
  return z
end

function nacb(x::narb, y::narb)
  z = nacb()
  ccall((:acb_set_arb_arb, libarb), Nothing,
        (Ref{nacb}, Ref{narb}, Ref{narb}),
        z, x, y)
  return z
end

function nacb(a::fmpz, p::Int)
  z = nacb()
  ccall((:acb_set_round_fmpz, libarb), Nothing,
        (Ref{nacb}, Ref{fmpz}, Int),
        z, a, p)
  return z  
end

function nacb_poly(c0::nacb, c1::Int)
  z = nacb_poly()
  ccall((:acb_poly_set_coeff_si, libarb), Nothing,
        (Ref{nacb_poly}, Int, Int),
        z, 1, c1)
  ccall((:acb_poly_set_coeff_acb, libarb), Nothing,
        (Ref{nacb_poly}, Int, Ref{nacb}),
        z, 0, c0)
  return z
end

function nacb_poly(c0::Int, c1::Int)
  z = nacb_poly()
  ccall((:acb_poly_set_coeff_si, libarb), Nothing,
        (Ref{nacb_poly}, Int, Int),
        z, 1, c1)
  ccall((:acb_poly_set_coeff_acb, libarb), Nothing,
        (Ref{nacb_poly}, Int, Ref{nacb}),
        z, 0, c0)
  return z
end

################################################################################

function native_string(x::narb)
  d = precision(x)
  d = clamp(d, 10, 200)
  d = trunc(Int, d*0.30103)
  cstr = ccall((:arb_get_str, Nemo.libarb), Ptr{UInt8},
         (Ref{narb}, Int, UInt),
         x, Int(d), UInt(0))
  r = unsafe_string(cstr)
  ccall((:flint_free, libflint), Nothing, (Ptr{UInt8},), cstr)
  return r
end

function Base.show(io::IO, a::narb)
  print(io, native_string(a))
end

function Base.show(io::IO, a::nacb)
  print(io, native_string(real(a)))
  print(io, " + i*")
  print(io, native_string(imag(a)))
end

function Base.show(io::IO, a::nacb_poly)
  first = true
  for i in 0:length(a)-1
    if !first
      print(io, " + ")
    end
    first = false
    print(io, "(")
    print(io, coeff(a, i))
    print(io, ")*ε^"*string(i))
  end
  if first
    print(io, "0")
  end
end

function Base.show(io::IO, a::nacb_mat)
  println(io, string(nrows(a))*" by "*string(ncols(a)))
  println(io, " [")
  for i in 1:nrows(a), j in 1:ncols(a)
    println(io, a[i,j])
  end
  println(io, "]")
end

#### boring arithmetic #########################################################

nrows(a::nacb_mat) = a.r

ncols(a::nacb_mat) = a.c

function getindex!(z::nacb, a::nacb_mat, i::Int, j::Int)
  GC.@preserve a begin
    aij = ccall((:acb_mat_entry_ptr, libarb), Ptr{nacb},
                (Ref{nacb_mat}, Int, Int),
                a, i - 1, j - 1)
    ccall((:acb_set, libarb), Nothing,
          (Ref{nacb}, Ptr{nacb}),
          z, aij)
  end
  return z
end

function Base.getindex(a::nacb_mat, i::Int, j::Int)
  @assert 0 < i <= nrows(a)
  @assert 0 < j <= ncols(a)
  z = nacb()
  getindex!(z, a, i, j)
  return z
end

function Base.setindex!(a::nacb_mat, i::Int, j::Int, b::nacb)
  GC.@preserve a begin
    aij = ccall((:acb_mat_entry_ptr, libarb), Ptr{nacb},
                (Ref{nacb_mat}, Int, Int),
                 a, i - 1, j - 1)
    ccall((:acb_set, libarb), Nothing,
          (Ptr{nacb}, Ref{nacb}),
          aij, b)
  end
  return b
end

function Base.length(a::nacb_poly)
  return a.length
end

function coeff(a::nacb_poly, n::Int)
  @assert n >= 0
  z = nacb()
  ccall((:acb_poly_get_coeff_acb, libarb), Nothing,
        (Ref{nacb}, Ref{nacb_poly}, Int),
        z, a, n)
  return z
end

function setcoeff!(a::nacb_poly, n::Int, b::nacb)
  @assert n >= 0
  ccall((:acb_poly_set_coeff_acb, libarb), Nothing,
        (Ref{nacb_poly}, Int, Ref{nacb}),
        a, n, b)
  return a
end


function precision(x::narb)
  return ccall((:arb_rel_accuracy_bits, libarb), Int,
               (Ref{narb}, ),
               x)
end

function precision(z::nacb)
  return ccall((:arb_rel_accuracy_bits, libarb), Int,
               (Ref{nacb}, ),
               z)
end

function precision_inc(a::Union{narb, nacb}, b::Int)
  p = precision(a)
  return clamp(p, 1, ARF_PREC_EXACT - b) + b
end

function Base.real(a::nacb)
  z = narb()
  ccall((:acb_get_real, libarb), Nothing,
        (Ref{narb}, Ref{nacb}),
        z, a)
  return z
end

function Base.imag(a::nacb)
  z = narb()
  ccall((:acb_get_imag, libarb), Nothing,
        (Ref{narb}, Ref{nacb}),
        z, a)
  return z
end


function convert(::Type{Float64}, x::narb)
  GC.@preserve x begin
    t = ccall((:arb_mid_ptr, libarb), Ptr{arf_struct}, (Ref{narb}, ), x)
    return ccall((:arf_get_d, libarb), Float64, (Ptr{arf_struct}, Int), t, 4)
  end
end

function convert(::Type{ComplexF64}, a::nacb)
  GC.@preserve a begin
    r = ccall((:acb_real_ptr, libarb), Ptr{narb}, (Ref{nacb}, ), a)
    i = ccall((:acb_imag_ptr, libarb), Ptr{narb}, (Ref{nacb}, ), a)
    t = ccall((:arb_mid_ptr, libarb), Ptr{arf_struct}, (Ptr{narb}, ), r)
    u = ccall((:arb_mid_ptr, libarb), Ptr{arf_struct}, (Ptr{narb}, ), i)
    v = ccall((:arf_get_d, libarb), Float64, (Ptr{arf_struct}, Int), t, 4)
    w = ccall((:arf_get_d, libarb), Float64, (Ptr{arf_struct}, Int), u, 4)
  end
  return complex(v, w)
end


function numerator!(z::fmpz, x::fmpq)
  ccall((:fmpq_numerator, libflint), Nothing,
        (Ref{fmpz}, Ref{fmpq}),
        z, x)
  return z
end

function denominator!(z::fmpz, x::fmpq)
  ccall((:fmpq_denominator, libflint), Nothing,
        (Ref{fmpz}, Ref{fmpq}),
        z, x)
  return z
end

function Base.zero(::Type{nacb})
  return nacb()
end

function zero!(z::fmpq)
  ccall((:fmpq_set_si, libflint), Nothing,
        (Ref{fmpq}, Int, UInt),
        z, Int(1), UInt(1))
  return z
end

function zero!(z::narb)
  ccall((:arb_zero, libarb), Nothing,
        (Ref{narb}, ),
        z)
  return z
end

function zero!(z::nacb)
  ccall((:acb_zero, libarb), Nothing,
        (Ref{nacb}, ),
        z)
  return z
end


function Base.one(::Type{nacb})
  return one!(nacb())
end

function Base.one(::Type{nacb_poly})
  z = nacb_poly()
  ccall((:acb_poly_set_si, libarb), Nothing,
        (Ref{nacb_poly}, Int),
        z, 1)
  return z  
end

function one!(z::fmpq)
  ccall((:fmpq_set_si, libflint), Nothing,
        (Ref{fmpq}, Int, UInt),
        z, Int(1), UInt(1))
  return z
end

function one!(z::narb)
  ccall((:arb_set_ui, libarb), Nothing,
        (Ref{narb}, UInt),
        z, UInt(1))
  return z
end

function one!(z::nacb)
  ccall((:acb_set_ui, libarb), Nothing,
        (Ref{nacb}, UInt),
        z, UInt(1))
  return z
end


function set!(z::narb, x::fmpq, p::Int)
  ccall((:arb_set_fmpq, libarb), Nothing,
        (Ref{narb}, Ref{fmpq}, Int),
        z, x, p)
  return z
end

function add!(z::fmpq, x::fmpq, y::Int)
  ccall((:fmpq_add_si, libflint), Nothing,
        (Ref{fmpq}, Ref{fmpq}, Int),
        z, x, y)
end

function mul!(z::fmpq, x::fmpq, y::fmpq)
  ccall((:fmpq_mul, libflint), Nothing,
        (Ref{fmpq}, Ref{fmpq}, Ref{fmpq}),
        z, x, y)
  return z
end

function mul!(z::narb, x::narb, y::Int, p::Int)
  ccall((:arb_mul_si, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int, Int),
        z, x, y, p)
  return z
end

function mul!(z::narb, x::narb, y::fmpz, p::Int)
  ccall((:arb_mul_fmpz, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{fmpz}, Int),
        z, x, y, p)
  return z
end

function mul!(z::nacb, x::nacb, y::fmpz, p::Int)
  ccall((:acb_mul_fmpz, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{fmpz}, Int),
        z, x, y, p)
  return z
end

function mul!(z::narb, x::narb, y::narb, p::Int)
  ccall((:arb_mul, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, x, y, p)
  return z
end

function mul!(z::nacb, x::nacb, y::nacb, p::Int)
  ccall((:acb_mul, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{nacb}, Int),
        z, x, y, p)
  return z
end

function div!(z::fmpq, x::fmpq, y::fmpq)
  ccall((:fmpq_div, libflint), Nothing,
        (Ref{fmpq}, Ref{fmpq}, Ref{fmpq}),
        z, x, y)
  return z
end

function div!(z::narb, x::narb, y::Int, p::Int)
  ccall((:arb_div_si, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int, Int),
        z, x, y, p)
  return z
end

function div!(z::narb, x::narb, y::fmpz, p::Int)
  ccall((:arb_div_fmpz, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{fmpz}, Int),
        z, x, y, p)
  return z
end

function div!(z::nacb, x::nacb, y::fmpz, p::Int)
  ccall((:acb_div_fmpz, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{fmpz}, Int),
        z, x, y, p)
  return z
end

function mul!(z::T, x::T, y::fmpq, p::Int) where T <: Union{narb, nacb}
  t = fmpz()
  numerator!(t, y)
  mul!(z, x, t, p)
  denominator!(t, y)
  div!(z, z, t, p)
  return z
end

function mullow(a::nacb_poly, b::nacb_poly, ord::Int, p::Int)
  z = nacb_poly()
  ccall((:acb_poly_mullow, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Ref{nacb_poly}, Int, Int),
        z, a, b, ord, p)
  return z  
end

function add!(z::narb, x::narb, y::Int, p::Int)
  ccall((:arb_add_si, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int, Int),
        z, x, y, p)
  return z
end

function add!(z::narb, x::narb, y::narb, p::Int)
  ccall((:arb_add_si, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, x, y, p)
  return z
end

function add!(z::nacb, x::nacb, y::Int, p::Int)
  ccall((:acb_add_si, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int, Int),
        z, x, y, p)
  return z
end

function add!(z::nacb, x::nacb, y::nacb, p::Int)
  ccall((:acb_add, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{nacb}, Int),
        z, x, y, p)
  return z
end

function sub!(z::fmpq, x::fmpq, y::fmpq)
  ccall((:fmpq_sub, libflint), Nothing,
        (Ref{fmpq}, Ref{fmpq}, Ref{fmpq}),
        z, x, y)
  return z
end

function sub!(z::narb, x::narb, y::Int, p::Int)
  ccall((:arb_sub_si, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int, Int),
        z, x, y, p)
  return z
end

function sub!(z::nacb, x::nacb, y::Int, p::Int)
  ccall((:acb_sub_si, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int, Int),
        z, x, y, p)
  return z
end

function sub!(z::nacb, x::nacb, y::nacb, p::Int)
  ccall((:acb_sub, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{nacb}, Int),
        z, x, y, p)
  return z
end

function neg!(z::nacb, a::nacb)
  ccall((:acb_neg, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}),
        z, a)
  return z
end

function neg!(z::nacb_poly, a::nacb)
  ccall((:acb_poly_neg, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}),
        z, a)
  return z
end

function neg!(z::Union{narb, nacb, nacb_poly})
    neg!(z, z)
end

function sub!(z::nacb, x::Int, y::nacb, p::Int)
  sub!(z, y, x, p)
  neg!(z, z)
  return z
end

function div!(z::nacb, x::nacb, y::nacb, p::Int)
  ccall((:acb_div, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{nacb}, Int),
        z, x, y, p)
  return z
end

function div!(z::nacb, x::nacb, y::Int, p::Int)
  ccall((:acb_div_si, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int, Int),
        z, x, y, p)
  return z
end

function Base.inv(x::nacb, p::Int)
  z = nacb()
  ccall((:acb_inv, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int),
        z, x, p)
  return z
end

function sqrt!(z::nacb, x::nacb, p::Int)
  ccall((:acb_sqrt, libarb), Nothing,
      (Ref{nacb}, Ref{nacb}, Int),
      z, x, p)
  return z
end

function log!(z::nacb, x::nacb, p::Int)
  ccall((:acb_log, libarb), Nothing,
      (Ref{nacb}, Ref{nacb}, Int),
      z, x, p)
  return z
end

function pow!(z::nacb, x::nacb, y::Int, p::Int)
  ccall((:acb_pow_si, libarb), Nothing,
      (Ref{nacb}, Ref{nacb}, Int, Int),
      z, x, y, p)
  return z
end

function pow!(z::nacb, x::nacb, y::fmpq, p::Int)
  Y = narb(y, p)
  ccall((:acb_pow_arb, libarb), Nothing,
      (Ref{nacb}, Ref{nacb}, Ref{narb}, Int),
      z, x, Y, p)
  return z
end


function ldexp!(z::nacb, x::nacb, y::Int)
  ccall((:acb_mul_2exp_si, libarb), Nothing,
      (Ref{nacb}, Ref{nacb}, Int),
      z, x, y)
  return z
end

function Base.:(+)(x::Int, y::nacb)
  return add!(nacb(), y, x, precision_inc(y, 60))
end

function Base.:(+)(y::nacb, x::Int)
  return add!(nacb(), y, x, precision_inc(y, 60))
end

function Base.:(-)(x::Int, y::nacb)
  return sub!(nacb(), x, y, precision_inc(y, 60))
end

function Base.:(-)(x::nacb, y::Int)
  return sub!(nacb(), x, y, precision_inc(x, 60))
end

function mul(x::nacb, y::Union{fmpq, nacb}, p::Int)
  return mul!(nacb(), x, y, p)
end

function Base.sqrt(x::nacb, p::Int)
  return sqrt!(nacb(), x, p)
end

function Base.log(x::nacb, p::Int)
  return log!(nacb(), x, p)
end

function Base.div(x::nacb, y::nacb, p::Int)
  return div!(nacb(), x, y, p)
end

function pow(x::nacb, y::Union{Int, fmpz, fmpq, narb, nacb}, p::Int)
  return pow!(nacb(), x, y, p)
end

function solve(a::nacb_mat, b::nacb_mat, p::Int)
  z = nacb_mat(ncols(a), ncols(b))
  ok = ccall((:acb_mat_solve, libarb), Cint,
             (Ref{nacb_mat}, Ref{nacb_mat}, Ref{nacb_mat}, Int),
             z, a, b, p)
  @assert ok != 0
  return z
end

function mul(a::nacb_mat, b::nacb_mat, p::Int)
  z = nacb_mat(nrows(a), ncols(b))
  ccall((:acb_mat_mul, libarb), Cint,
        (Ref{nacb_mat}, Ref{nacb_mat}, Ref{nacb_mat}, Int),
        z, a, b, p)
  return z
end

function Base.hcat(a::nacb_mat, b::nacb_mat)
  @assert nrows(a) == nrows(b)
  z = nacb_mat(nrows(a), ncols(a) + ncols(b))
  t = nacb()
  for i in 1:nrows(a)
    for j in 1:ncols(a)
      getindex!(t, a, i, j)
      setindex!(z, i, j, t)
    end
    for j in 1:ncols(b)
      getindex!(t, b, i, j)
      setindex!(z, i, ncols(a) + j, t)
    end
  end
  return z
end

function overlaps(a::nacb, b::nacb)
  r = ccall((:acb_overlaps, libarb), Cint,
            (Ref{nacb}, Ref{nacb}),
            a, b)
  return r != 0
end

function swap!(a::nacb, b::nacb)
  ccall((:acb_swap, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}),
        a, b)
end



#### special functions required for field F with elem_type T ###################

# integer check
function isinteger_with_integer(x::fmpq)
  if isone(denominator(x))
    return true, numerator(x)
  else
    return false, ZZ()
  end
end

# ball approximation
function narb(x::fmpq, p::Int)
  z = narb()
  ccall((:arb_set_fmpq, libarb), Nothing,
      (Ref{narb}, Ref{fmpq}, Int),
      z, x, p)
  return z
end

function nacb(x::fmpq, p::Int)
  z = nacb()
  ccall((:acb_set_fmpq, libarb), Nothing,
      (Ref{nacb}, Ref{fmpq}, Int),
      z, x, p)
  return z
end


#### diff equation satisfied by pFq ############################################

# x is the independent variable. θ = x d/dx is the operator.
# θ*x = x*θ - x

# hack
function Nemo.AbstractAlgebra.Generic.derivative(f::FracElem)
  n = numerator(f)
  d = denominator(f)
  return (derivative(n)*d - n*derivative(d))//d^2
end

# equation c maintains θ on right as c0(x) + c1(x)*θ + c2(x)*θ^2 + ...
# do c = fmul*θ*fmul^-1*c + d*c where logdm = fmul'/fmul
function equ_θ_mul_add!(c::Vector, logdm, arg, d)
  FFx = parent(arg)   # F(x)
  x = FFx(gen(base_ring(FFx)))
  n = length(c)
  push!(c, zero(FFx))
  v = divexact(arg, derivative(arg))
  for i in n:-1:0
    t = (i == 0) ? FFx(0) : divexact(c[1 + i - 1], x)
    t += derivative(c[1 + i]) - c[1 + i]*logdm
    c[1 + i] = c[1 + i]*d + v*t
  end
end

# equation c maintains θ on left as c0(x) + θ*c1(x) + θ^2*c2(x) + ...
# do c = c*θ + d
function equ_mul_θ_add!(c::Vector, d, Fx)
  n = length(c)
  push!(c, Fx())
  for i in n:-1:0
    c[1 + i] = (i == n) ? Fx(0) : -shift_left(derivative(c[1 + i]), 1)
    c[1 + i] += (i == 0) ? d : c[1 + i - 1]
  end
end

# return either c1 - c2 or c2 - c1, doesn't matter which
function equ_sub!(c1::Vector, c2::Vector)
  if length(c1) < length(c2)
    (c1, c2) = (c2, c1)
  end
  for i in 1:length(c2)
    c1[i] -= c2[i]
  end
  return c1
end

# Return equation satisfied by m(x)*pFq(arg(x)) with θ on the left.
# Allowed m's are those for which logdm(x) = m'(x)/m(x) is rational.
# Equation is returned as the dense coefficient list of powers of x.
function hyp_equ(a::Vector{T}, b::Vector{T}, logdm, arg) where T
  FFx = parent(arg)   # F(x)
  Fx = base_ring(FFx) # F[x]
  F = base_ring(Fx)   # F
  @assert FFx == parent(logdm)

  # c = equation with θ on the right
  c = elem_type(FFx)[arg]
  for ai in a
    equ_θ_mul_add!(c, logdm, arg, ai - 1)
  end
  c2 = elem_type(FFx)[one(FFx)]  
  equ_θ_mul_add!(c2, logdm, arg, F(0)) # "extra" denominator param 1
  for bi in b
    equ_θ_mul_add!(c2, logdm, arg, bi - 1)
  end
  c = equ_sub!(c, c2)

  # cancel content while θ is on the right
  c = map(p -> (numerator(p), denominator(p)), c)
  cont = zero(Fx)
  den = one(Fx)
  for ci in c
    cont = gcd(cont, ci[1])
    den = lcm(den, ci[2])
  end
  c = elem_type(Fx)[divexact(ci[1], cont)*divexact(den, ci[2]) for ci in c]

  # move θ to the left
  tP = elem_type(Fx)[]
  for ci in reverse(c)
    equ_mul_θ_add!(tP, ci, Fx)
  end

  # transpose so θ is inside P0(θ) + P1(θ)*x + ... + Pn(θ)*x^n
  Fθ, θ = PolynomialRing(F, "θ")
  P = elem_type(Fθ)[]
  for i in (length(tP) - 1):-1:0
    for j in 0:degree(tP[1 + i])
      while j >= length(P)
        push!(P, zero(Fθ))
      end
      setcoeff!(P[1 + j], i, coeff(tP[1 + i], j))
    end
  end

  return P
end

#### evaluation of series solutions to diff equations ##########################

function get_multiplicity!(ns::Vector{fmpz}, iv::Vector{Vector{T}}, n::Int) where T
  @assert length(ns) == length(iv)
  nextu = Vector{T}[]
  while !isempty(ns) && ns[1] == n
    popfirst!(ns)
    push!(nextu, popfirst!(iv))
  end
  return nextu
end

# Return δ by ν matrix M such that for any vector iv of ν initial values the
# solution f(z) determined by these initial values and its derivatives is M.iv.
# Currently doesn't work well when δ > 1 and z contains 0. This case is
# fortunately not needed.
function eval_basis(
  P::Vector,
  λ::T, ns::Vector{fmpz}, # indicial roots are λ+ns[1], λ+ns[2], ..., λ+ns[ν] 
  z::nacb,
  δ::Int,       # calculate f(z), f'(z), ... f^(δ-1)(z)
  normal::Bool, # continuity along z < 0 is determined by counterclockwise direction
  invz::nacb,   # 1/z only used if normal = false
  p::Int) where T

  ν = length(ns)    # number of initial conditions
  τ = 0             # strict bound on the power of log(z) thus far
  maxN = 10 + 20*p  # max number of terms to sum

  @assert length(ns) > 0
  @assert ns[1] == 0

  s = length(P) - 1
  @assert s > 0

  Fθ = parent(P[1])
  θ = gen(Fθ)
  F = base_ring(Fθ)

  # u is an array of the last s coefficients
  # each solution is an array of coefficients of log^k(z)/k!  0 <= k < τ
  # each coefficient is an array of ν elements of F
  u = [Vector{T}[] for i in 1:s]

  # Σ is an array of terms for the i^th derivative
  # each terms is an array of coefficients of log^k(z)/k!  0 <= k < τ
  # each coefficient is an array of ν elements of nacb
  Σ = [Vector{nacb}[] for i in 1:δ]

  #
  ns = deepcopy(ns)
  iv = [[i == j ? one(F) : zero(F) for i in 1:ν] for j in 1:ν]

  t = nacb()
  zn = one(nacb)

  changedcount = 0
  n = 0
  while n < maxN
    # evaluate the next coefficient u_n
    Pn = map(a -> a(θ + (λ + n)), P)
    un = get_multiplicity!(ns, iv, n)
    rhs = [[zero(F) for k in 1:ν] for j in 1:τ]
    for i in 1:s
      # Pn[1 + i] is a polynomial in θ. Apply it to u[i] where θ is
      # the shift operator and sub result from rhs
      for l in 0:degree(Pn[1 + i])
        for j in 1:length(u[i])-l, k in 1:ν
          sub!(rhs[j][k], rhs[j][k], coeff(Pn[1+i], l)*u[i][j + l][k])
        end
      end
    end
    # μ = number of initial conditions at λ + n = multiplicity of root
    μ = length(un)
    for i in 0:μ-1
      @assert iszero(coeff(Pn[1], i))
    end
    for j in 1:τ
      push!(un, rhs[j])
    end
    for i in 1:τ
      for j in 1:i-1
        for k in 1:ν
          un[1+μ+τ-i][k] -= coeff(Pn[1], μ+j)*un[1+μ+τ-i+j][k]
        end
      end
      for k in 1:ν
        un[1+μ+τ-i][k] //= coeff(Pn[1], μ)
      end
    end
    # trim zeros off un
    while length(un) > τ && all(iszero, un[end])
      pop!(un)
    end
    τ = max(τ, length(un))
    pop!(u)
    pushfirst!(u, un)

    # add un*z^n to sum
    changed = false
    znrow = [[mul(zn, un[1+i][k], p) for k in 1:ν] for i in 0:τ-1]
    for d in 0:δ-1
      if d > 0
        # differentiate un
        for i in 0:τ-1, k in 1:ν          
          mul!(znrow[1+i][k], znrow[1+i][k], λ+n-d+1, p)
          if i+1 < τ
            add!(znrow[1+i][k], znrow[1+i][k], znrow[1+i+1][k], p)
          end
        end
      end
      while length(Σ[1+d]) < τ
        push!(Σ[1+d], [nacb() for k in 1:ν])
      end
      for i in 0:τ-1, k in 1:ν
        add!(t, Σ[1+d][1+i][k], znrow[1+i][k], p)
        changed = changed || !overlaps(Σ[1+d][1+i][k], t)
        swap!(Σ[1+d][1+i][k], t)
      end
    end

    mul!(zn, zn, z, p)
    n += 1
    if !changed
      if (changedcount += 1) > 1
        break
      end
    else
        changedcount = 0
    end
  end

  # TODO add a tail bound

  # evaluate the polynomials in log(z)
  logz = normal ? log(z, p) : neg!(log(invz, p))
  zλ = normal ? pow(z, λ, p) : pow(invz, -λ, p)
  M = nacb_mat(δ, ν)
  for d in 0:δ-1, k in 1:ν
    i = τ-1
    if i < 0
      t = zero(nacb)
    else
      t = Σ[1+d][1+i][k]
      while (i -= 1) >= 0
        div!(t, t, i + 1, p)
        mul!(t, t, logz, p)
        add!(t, t, Σ[1+d][1+i][k], p)
      end
    end
    mul!(t, t, d == 0 ? zλ : mul(zλ, pow(z, -d, p), p), p)
    setindex!(M, 1+d, k, t)
  end
  return M
end


function eval_bases(P, ρ, z, n, normal, invz, p)
  M = eval_basis(P, ρ[1][1], ρ[1][2], z, n, normal, invz, p)
  for i in 2:length(ρ)
    M = hcat(M, eval_basis(P, ρ[i][1], ρ[i][2], z, n, normal, invz, p))
  end
  return M
end


#### nFn-1 evaluation with parameters in a field F #############################

# return array of (λ_i, [n_i1 = 0, n_i2, ...], Ind_i) such that
#   the λ_i are distinct mod 1
#   0 = n_i1 <= n_i2 <= n_i3 .... <= n_ik are integers
#   for each i, the unordered list {λ_i + n_ij}_j = the unordered list {a[Ii[j]]}_j
#   the Ind_i[j] are of course all distinct
function partition_mod1(a::Vector{T}) where T
  F = parent(a[1])
  r = Tuple{elem_type(F), Vector{elem_type(ZZ)}, Vector{Int}}[]
  for i in 1:length(a)
    found = false
    for j in 1:length(r)
      ok, k = isinteger_with_integer(a[i] - r[j][1])
      if ok
        push!(r[j][2], k)
        push!(r[j][3], i)
        found = true
        break
      end
    end
    if !found
      push!(r, (a[i], [ZZ(0)], [i]))
    end
  end
  for j in 1:length(r)
    sort!(r[j][2])
    r[j] = (r[j][1] + r[j][2][1], r[j][2] .- r[j][2][1], r[j][3]) 
  end
  return r
end

# series in a + b*ε
function rgamma_series(a::nacb, b::Int, ord::Int, p::Int)
  z = nacb_poly(a, b)
  ccall((:acb_poly_rgamma_series, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Int, Int),
        z, z, ord, p)
  return z
end

function gamma_series(a::nacb, b::Int, ord::Int, p::Int)
  z = nacb_poly(a, b)
  ccall((:acb_poly_gamma_series, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Int, Int),
        z, z, ord, p)
  return z
end

function rising_factorial_series(a::nacb, b::Int, t::fmpz, ord::Int, p::Int)
  @assert t >= 0
  z = nacb_poly(a, b)
  # TODO handle large t
  ccall((:acb_poly_rising_ui_series, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, UInt, Int, Int),
        z, z, UInt(t), ord, p)
  return z
end

# π*ε*csc(x*π*ε)
function cscpi_series(x::fmpz, ord::Int, p::Int)
  z = nacb_poly()
  setcoeff!(z, 0, nacb(1//x, p))
  t = narb()
  for n in 1:div(ord-1,2)
    ccall((:arb_const_pi, libarb), Nothing,
          (Ref{narb}, Int),
          t, p)
    pow!(t, t, 2*n, p)
    mul!(t, t, bernoulli(2*n)/factorial(ZZ(2*N))*
               (-1)^(n+1)*(4^n-2)*x^(2*n-1), p)
    setcoeff!(z, 2*n, nacb(t))
  end
  return z
end

function gamma(a::nacb, p::Int)
  z = nacb()
  ccall((:acb_gamma, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int),
        z, a, p)
  return z
end

function rgamma(a::nacb, p::Int)
  z = nacb()
  ccall((:acb_rgamma, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int),
        z, a, p)
  return z
end

# the numerator parameters are λ + am[1], ..., λ + am[m], as[1], ... , as[.]
# produce the initial values for the exponents λ + am[1], ..., λ + am[m]
# computed by replacing λ + am[1+i] -> λ + am[1+i] + i*ε and letting ε -> 0
# no difference λ - as[i] is an integer
function hyp_initial_values_at_∞(
  λ::T, am::Vector{fmpz}, as::Vector{T},
  b::Vector{T},
  p::Int) where T

  m = length(am)
  r = [zero(nacb) for j in 1:m]
  j = 0
  while j < m
    n = am[1 + j]
    # μ = multiplicity of exponent λ + n
    μ = 1
    while j + μ < m && am[1 + j + μ] == n
      μ += 1
    end
    # compute coeffs r[1+j], ..., r[1+j+μ-1]
    # of z^(λ+n), ..., z^(λ+n)*log^(μ-1)(z)/(μ-1)!
    for k in 0:m-1
      # t is the index of the term in the pFq series
      t = n - am[1 + k]
      if t < 0
        continue
      end
      # compute the order required for the power series O(ε^ord)
      ord = 1
      for i in 0:m-1
        if i != k && am[1 + i] - am[1 + k] - t <= 0
          # a gamma factor has a simple pole
          ord += 1
        end
      end
      s = one(nacb_poly)
      # handle integer differences between numerator parameters
      for i in 0:m-1
        if i == k
          s2 = rising_factorial_series(nacb(λ + am[1 + i], p), i, t, ord, p)
          s = mullow(s, s2, ord, p)
        else
          s2 = rgamma_series(nacb(λ + am[1 + i], p), i, ord, p)
          s = mullow(s, s2, ord, p)
          d = am[1 + i] - am[1 + k] - t
          if d > 0
            s2 = gamma_series(nacb(d, p), i - k, ord, p)
            s = mullow(s, s2, ord, p)
          else
            s2 = rgamma_series(nacb(1 - d, p), k - i, ord, p)
            s2 = mullow(s2, cscpi_series(ZZ(i - k), ord, p), ord, p)
            s = mullow(s, s2, ord, p)
            isodd(d) && neg!(s)
          end
        end
      end
      # no difference A - λ - am[1 + i] is an integer
      for A in as
        s2 = gamma_series(nacb(A - λ - am[1 + k] - t, p), -k, ord, p)
        s = mullow(s, s2, ord, p)
      end
      # the difference B - λ - am[1 + i] could be an integer, but it's rgamma
      for B in b
        s2 = rgamma_series(nacb(B - λ - am[1 + k] - t, p), -k, ord, p)
        s = mullow(s, s2, ord, p)
      end
      # now determine r[1+j], ..., r[1+j+μ-1]. They are already zeroed.
      for i in 0:min(ord, μ) - 1
        f = 1//factorial(t) # TODO
        if i > 0
          f *= ZZ(k)^i//factorial(ZZ(i))
        end
        ff = mul(coeff(s, ord - 1 - i), f, p)
        if iseven(t)
          add!(r[1 + j + i], r[1 + j + i], ff, p)
        else
          sub!(r[1 + j + i], r[1 + j + i], ff, p)
        end
      end
    end
    j += μ
  end
  f = one(nacb)
  for A in as
    mul!(f, f, rgamma(nacb(A, p), p), p)
  end
  for B in b
    mul!(f, f, gamma(nacb(B, p), p), p)
  end
  iv = nacb_mat(m, 1)
  for i in 1:m
    setindex!(iv, i, 1, mul(r[i], f, p))
  end
  return iv
end

function compute_f_at_0(a::Vector{T}, b::Vector{T}, z::nacb, p::Int) where T
  F = parent(a[1])
  Fx, x = PolynomialRing(F, "x")
  P = hyp_equ(a, b, Fx(0)//Fx(1), x//Fx(1))
  e = eval_basis(P, F(0), [ZZ(0)], z, 1, true, nacb(), p)
  return e[1, 1]  
end

function compute_f_at_∞(a::Vector{T}, b::Vector{T}, z::nacb, p::Int) where T
  F = parent(a[1])
  Fx, x = PolynomialRing(F, "x")
  # ∞ -> z
  P = hyp_equ(a, b, Fx(0)//Fx(1), -Fx(1)//x)
  nz = neg!(nacb(), z)
  rz = inv(nz, p)
  s = zero(nacb)
  for (λ, ns, indices) in partition_mod1(a)
    arest = T[a[j] for j in 1:length(a) if !(j in indices)]
    iv = hyp_initial_values_at_∞(λ, ns, arest, b, p)
    e = eval_basis(P, λ, ns, rz, 1, false, nz, p)
    add!(s, s, mul(e, iv, p)[1, 1], p)
  end
  return s
end

function compute_f_at_1(a::Vector{T}, b::Vector{T}, z::nacb, p::Int) where T
  n = length(a)
  F = parent(a[1])
  Fx, x = PolynomialRing(F, "x")
  w = nacb(QQ(3//4), p)
  # 0 -> w
  P = hyp_equ(a, b, Fx(0)//Fx(1), x//Fx(1))
  e0 = eval_basis(P, F(0), [ZZ(0)], w, n, true, nacb(), p)
  for i in 2:2:n
    setindex!(e0, i, 1, neg!(e0[i, 1])) # correct odd derivatives
  end
  # roots of indicial equation are 0:n-2 and σ
  σ::T = sum(b) - sum(a)
  ok, ZZσ = isinteger_with_integer(σ)
  if n < 2
    ρ = [(σ, [ZZ(0)])]
  elseif ok
    roots = sort!(push!(ZZ.(0:n-2), ZZσ))
    ρ = [(F(roots[1]), roots .- roots[1])]
  else
    ρ = [(σ, [ZZ(0)]), (F(0), ZZ.(0:n-2))]
  end
  # 1 -> w and 1 -> z
  P = hyp_equ(a, b, Fx(0)//Fx(1), (1 - x)//Fx(1))
  e1 = eval_bases(P, ρ, 1 - z, 1, true, nacb(), p)
  e2 = eval_bases(P, ρ, 1 - w, n, true, nacb(), p)
  # 0 -> w -> 1 -> z
  return mul(e1, solve(e2, e0, p), p)[1, 1]
end

function compute_f_anywhere(a::Vector{T}, b::Vector{T}, z::nacb, p::Int) where T
  F = parent(a[1])
  Fx, x = PolynomialRing(F, "x")
  # 0 -> (1-sqrt(1-z))/(1+sqrt(1-z))
  P = hyp_equ(a, b, -2*Fx(a[1])//(1 + x), 4*x//(1 + x)^2)
  t = sqrt(1 - z, p)
  s = 1 + t
  e = eval_basis(P, F(0), [ZZ(0)], div(1 - t, s, p), 1, true, nacb(), p)
  ldexp!(s, s, -1)
  return mul(pow(s, -2*a[1], p), e[1, 1], p)
end

# The non-regularized function is undefined if any bi is in {0, -1, -2, ...}.
# The regularized function is not implemented and may be defined via the
# initial values at the exponents 0 and the 1 - bi:
#   - any non-integer exponent has value 0
#   - only the greatest integer exponent n has nonzero value and this value is
#      prod_i (ai)_n / prod_i Γ(bi + n) * z^n / n!  (with no log(z))
function compute_pfq(a::Vector{T}, b::Vector{T}, z::nacb, p::Int) where T
  @assert length(a) == length(b) + 1
  zz = convert(Complex{Float64}, z)
  if abs2(zz) < 0.8
    return compute_f_at_0(a, b, z, p)
  elseif abs2(zz) > 1.2
    return compute_f_at_∞(a, b, z, p)
  elseif abs2(zz) < 2*real(zz) - 0.75
    return compute_f_at_1(a, b, z, p)
  else
    return compute_f_anywhere(a, b, z, p)  
  end
end

#### tests #####################################################################

@show compute_pfq([QQ(1//3), QQ(2//3)],
                  [QQ(2)],
                  nacb(QQ(-21//20), 100),
                  100)

@show compute_pfq([QQ(1//3), QQ(2//3)],
                  [QQ(2)],
                  nacb(QQ(21//20), 100),
                  100)

@show compute_pfq([QQ(1//3), QQ(2//3)],
                  [QQ(2)],
                  nacb(QQ(19//20), 100),
                  100)

@show compute_pfq([QQ(1//3), QQ(2//3)],
                  [QQ(2)],
                  nacb(QQ(-31//20), 100),
                  100)

@show compute_pfq([QQ(1//3), QQ(2//3)],
                  [QQ(2)],
                  nacb(QQ(31//20), 100),
                  100)


@show compute_pfq([QQ(1//3), QQ(2//3), QQ(1//2)],
                  [QQ(2), QQ(4//5)],
                  nacb(QQ(-21//20), 100),
                  100)

@show compute_pfq([QQ(1//3), QQ(2//3), QQ(1//2)],
                  [QQ(2), QQ(4//5)],
                  nacb(QQ(21//20), 100),
                  100)

@show compute_pfq([QQ(1//3), QQ(2//3), QQ(1//2)],
                  [QQ(2), QQ(4//5)],
                  nacb(QQ(19//20), 100),
                  100)

@show compute_pfq([QQ(1//3), QQ(2//3), QQ(1//2)],
                  [QQ(2), QQ(4//5)],
                  nacb(QQ(-31//20), 100),
                  100)

@show compute_pfq([QQ(1//3), QQ(2//3), QQ(1//2)],
                  [QQ(2), QQ(4//5)],
                  nacb(QQ(31//20), 100),
                  100)

nothing
