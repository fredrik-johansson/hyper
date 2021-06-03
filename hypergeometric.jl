using Nemo
import Nemo: libarb, libflint, arf_struct, arb_struct, expressify
import Nemo: zero!, one!, add!, sub!, mul!, div!, inv!, coeff,
       setcoeff!, mullow, gamma, rgamma, solve, overlaps, sub

const ARF_PREC_EXACT = typemax(Int)

#### five types: arb, acb, arb_poly, acb_poly, and acb_mat #####################

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


mutable struct narb_poly
  coeffs::Ptr{Nothing}
  length::Int
  alloc::Int

  function narb_poly()
    z = new()
    ccall((:arb_poly_init, libarb), Nothing, (Ref{narb_poly}, ), z)
    finalizer(_arb_poly_clear_fn, z)
    return z
  end
end

function _arb_poly_clear_fn(x::narb_poly)
  ccall((:arb_poly_clear, libarb), Nothing, (Ref{narb_poly}, ), x)
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

#### constructors ##############################################################

function narb(a::Int)
  z = narb()
  ccall((:arb_set_si, libarb), Nothing,
        (Ref{narb}, Int),
        z, a)
  return z
end

function nacb(a::Int)
  z = nacb()
  ccall((:acb_set_si, libarb), Nothing,
        (Ref{nacb}, Int),
        z, a)
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

function narb_poly(c0::narb)
  z = narb_poly()
  ccall((:arb_poly_set_coeff_arb, libarb), Nothing,
        (Ref{narb_poly}, Int, Ref{narb}),
        z, 0, c0)
  return z
end

function narb_poly(c0::narb, c1::Int)
  z = narb_poly()
  ccall((:arb_poly_set_coeff_si, libarb), Nothing,
        (Ref{narb_poly}, Int, Int),
        z, 1, c1)
  ccall((:arb_poly_set_coeff_arb, libarb), Nothing,
        (Ref{narb_poly}, Int, Ref{narb}),
        z, 0, c0)
  return z
end

function nacb_poly(c0::nacb, c1::nacb)
  z = nacb_poly()
  ccall((:acb_poly_set_coeff_acb, libarb), Nothing,
        (Ref{nacb_poly}, Int, Ref{nacb}),
        z, 1, c1)
  ccall((:acb_poly_set_coeff_acb, libarb), Nothing,
        (Ref{nacb_poly}, Int, Ref{nacb}),
        z, 0, c0)
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

#### output ####################################################################

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

function Base.show(io::IO, a::Union{narb_poly, nacb_poly})
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

Nemo.AbstractAlgebra.needs_parentheses(x::Union{narb, nacb, nacb_poly, nacb_mat}) = true


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

function Base.length(a::Union{narb_poly, nacb_poly})
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

function coeff(a::narb_poly, n::Int)
  @assert n >= 0
  z = narb()
  ccall((:arb_poly_get_coeff_arb, libarb), Nothing,
        (Ref{narb}, Ref{narb_poly}, Int),
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

function setcoeff!(a::narb_poly, n::Int, b::narb)
  @assert n >= 0
  ccall((:arb_poly_set_coeff_arb, libarb), Nothing,
        (Ref{narb_poly}, Int, Ref{narb}),
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

function Base.zero(::Type{narb})
  return narb()
end

function Base.zero(::Type{nacb_poly})
  return nacb_poly()
end

function Base.zero(::Type{narb_poly})
  return narb_poly()
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


function Base.one(::Type{narb})
  return one!(narb())
end

function Base.one(::Type{nacb})
  return one!(nacb())
end

function Base.one(::Type{narb_poly})
  z = narb_poly()
  ccall((:arb_poly_set_si, libarb), Nothing,
        (Ref{narb_poly}, Int),
        z, 1)
  return z  
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

function Base.max(x::narb, y::narb, p::Int)
  z = narb()
  ccall((:arb_max, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, x, y, p)
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

function mul!(z::nacb_poly, x::nacb_poly, y::nacb, p::Int)
  ccall((:acb_poly_scalar_mul, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Ref{nacb}, Int),
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

function mullow(a::narb_poly, b::narb_poly, ord::Int, p::Int)
  z = narb_poly()
  ccall((:arb_poly_mullow, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Ref{narb_poly}, Int, Int),
        z, a, b, ord, p)
  return z  
end

function mullow(a::nacb_poly, b::nacb_poly, ord::Int, p::Int)
  z = nacb_poly()
  ccall((:acb_poly_mullow, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Ref{nacb_poly}, Int, Int),
        z, a, b, ord, p)
  return z  
end

function div_series(a::nacb_poly, b::nacb_poly, ord::Int, p::Int)
  z = nacb_poly()
  ccall((:acb_poly_div_series, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Ref{nacb_poly}, Int, Int),
        z, a, b, ord, p)
  return z  
end

function add!(z::nacb_poly, a::nacb_poly, b::nacb_poly, p::Int)
  ccall((:acb_poly_add, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}, Ref{nacb_poly}, Int),
        z, a, b, p)
  return z  
end

function add!(z::narb, x::narb, y::Int, p::Int)
  ccall((:arb_add_si, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int, Int),
        z, x, y, p)
  return z
end

function max!(z::narb, x::narb, y::narb, p::Int)
  ccall((:arb_max, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, x, y, p)
  return z
end

function add!(z::narb, x::narb, y::narb, p::Int)
  ccall((:arb_add, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, x, y, p)
  return z
end

function add(a::narb, b::narb, p::Int)
  return add!(narb(), a, b, p)
end

function sub!(z::narb, x::narb, y::narb, p::Int)
  ccall((:arb_sub, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Ref{narb}, Int),
        z, x, y, p)
  return z
end

function sub(a::narb, b::narb, p::Int)
  return sub!(narb(), a, b, p)
end

function add!(z::narb_poly, x::narb_poly, y::narb_poly, p::Int)
  ccall((:arb_poly_add, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Ref{narb_poly}, Int),
        z, x, y, p)
  return z
end

function add(a::narb_poly, b::narb_poly, p::Int)
  return add!(narb_poly(), a, b, p)
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

function sub!(z::nacb, a::nacb, b::nacb, p::Int)
  ccall((:acb_sub, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Ref{nacb}, Int),
        z, a, b, p)
  return z
end

function add!(z::narb_poly, a::narb_poly, b::Int, p::Int)
  ccall((:arb_poly_add_si, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Int, Int),
        z, a, b, p)
  return z
end

function add(a::narb_poly, b::Int, p::Int)
  return add!(narb_poly(), a, b, p)
end

function add(b::Int, a::narb_poly, p::Int)
  return add!(narb_poly(), a, b, p)
end


function sub!(z::narb_poly, a::Int, b::narb_poly, p::Int)
  return neg!(add!(z, b, -a, p)) # TODO
end

function sub(a::Int, b::narb_poly, p::Int)
  return sub!(narb_poly(), a, b, p)
end

function neg!(z::nacb, a::nacb)
  ccall((:acb_neg, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}),
        z, a)
  return z
end

function neg!(z::narb_poly, a::narb_poly)
  ccall((:arb_poly_neg, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}),
        z, a)
  return z
end

function neg!(z::nacb_poly, a::nacb_poly)
  ccall((:acb_poly_neg, libarb), Nothing,
        (Ref{nacb_poly}, Ref{nacb_poly}),
        z, a)
  return z
end

function neg!(z::Union{narb, nacb, narb_poly, acb_poly})
  return neg!(z, z)
end

function Base.:-(a::T) where T <: Union{narb, nacb, narb_poly, acb_poly}
  return neg!(T(), a)
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

function div!(z::narb, x::narb, y::Int, p::Int)
  ccall((:arb_div_si, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int, Int),
        z, x, y, p)
  return z
end

function div!(z::nacb, x::nacb, y::Int, p::Int)
  ccall((:acb_div_si, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int, Int),
        z, x, y, p)
  return z
end

function inv!(z::narb, x::narb, p::Int)
  ccall((:arb_inv, libarb), Nothing,
        (Ref{narb}, Ref{narb}, Int),
        z, x, p)
  return z
end

function inv!(z::nacb, x::nacb, p::Int)
  ccall((:acb_inv, libarb), Nothing,
        (Ref{nacb}, Ref{nacb}, Int),
        z, x, p)
  return z
end

function Base.inv(x::nacb, p::Int)
  return inv!(nacb(), x, p)
end

function Base.abs(x::nacb, p::Int)
  z = narb()
  ccall((:acb_abs, libarb), Nothing,
        (Ref{narb}, Ref{nacb}, Int),
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

function mul(x::narb, y::Union{fmpq, narb}, p::Int)
  return mul!(narb(), x, y, p)
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

function Base.div(x::narb, y::Int, p::Int)
  return div!(narb(), x, y, p)
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

function inv_series(a::narb_poly, ord::Int, p::Int)
  z = narb_poly()
  ccall((:arb_poly_inv_series, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Int, Int),
        z, a, ord, p)
  return z  
end

function log_series(a::narb_poly, ord::Int, p::Int)
  z = narb_poly()
  ccall((:arb_poly_log_series, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Int, Int),
        z, a, ord, p)
  return z  
end

function log1p_series(a::narb_poly, ord::Int, p::Int)
  z = narb_poly()
  ccall((:arb_poly_log1p_series, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Int, Int),
        z, a, ord, p)
  return z  
end

function exp_series(a::narb_poly, ord::Int, p::Int)
  z = narb_poly()
  ccall((:arb_poly_exp_series, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Int, Int),
        z, a, ord, p)
  return z
end

function pow_series(a::narb_poly, b::narb, ord::Int, p::Int)
  z = narb_poly()
  ccall((:arb_poly_pow_arb_series, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Ref{narb}, Int, Int),
        z, a, b, ord, p)
  return z  
end

function mul!(z::narb_poly, a::narb_poly, b::narb, p::Int)
  ccall((:arb_poly_scalar_mul, libarb), Nothing,
        (Ref{narb_poly}, Ref{narb_poly}, Ref{narb}, Int),
        z, a, b, p)
  return z  
end

function mul(a::narb_poly, b::narb, p::Int)
  return mul!(narb_poly(), a, b, p)
end


function Base.factorial(a::Int, p::Int)
  z = narb()
  ccall((:arb_fac_ui, libarb), Nothing,
        (Ref{narb}, UInt, Int),
        z, UInt(a), p)
  return z
end

function add_error!(z::nacb, a::narb)
  GC.@preserve z begin
    t = ccall((:acb_real_ptr, libarb), Ptr{narb}, (Ref{nacb}, ), z)
    ccall((:arb_add_error, libarb), Nothing, (Ptr{narb}, Ref{narb}), t, a)
    t = ccall((:acb_imag_ptr, libarb), Ptr{narb}, (Ref{nacb}, ), z)
    ccall((:arb_add_error, libarb), Nothing, (Ptr{narb}, Ref{narb}), t, a)
  end
  return z
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
    t = (i == 0) ? FFx(0) : divexact(c[1+i-1], x)
    t += derivative(c[1+i]) - c[1 + i]*logdm
    c[1+i] = c[1+i]*d + v*t
  end
end

# equation c maintains θ on left as c0(x) + θ*c1(x) + θ^2*c2(x) + ...
# do c = c*θ + d
function equ_mul_θ_add!(c::Vector, d, Fx)
  n = length(c)
  push!(c, Fx())
  for i in n:-1:0
    c[1+i] = (i == n) ? Fx(0) : -shift_left(derivative(c[1+i]), 1)
    c[1+i] += (i == 0) ? d : c[1+i-1]
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
    equ_θ_mul_add!(c, logdm, arg, ai-1)
  end
  c2 = elem_type(FFx)[one(FFx)]  
  equ_θ_mul_add!(c2, logdm, arg, F(0)) # "extra" denominator param 1
  for bi in b
    equ_θ_mul_add!(c2, logdm, arg, bi-1)
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
  Px = elem_type(Fx)[]
  for ci in reverse(c)
    equ_mul_θ_add!(Px, ci, Fx)
  end

  # transpose so θ is inside P0(θ) + P1(θ)*x + ... + Pn(θ)*x^n
  Fθ, θ = PolynomialRing(F, "θ")
  Pθ = elem_type(Fθ)[]
  for i in (length(Px)-1):-1:0
    for j in 0:degree(Px[1+i])
      while j >= length(Pθ)
        push!(Pθ, zero(Fθ))
      end
      setcoeff!(Pθ[1+j], i, coeff(Px[1+i], j))
    end
  end

  return (Pθ, Px)
end


#### evaluation of series solutions to diff equations ##########################
#=
The setup follows Mezzarobba 2019 Truncation Bounds

With θ on the left we have a differential operator

  P(x,θ) = θ^r*Px_r + ... + θ*Px_1 + Px_0,  Px_i in F[x]
         = Pθ_s*x^s + ... + Pθ_1*x + Pθ_0,  Pθ_i in F[θ]

The radius of convergence of the solutions about 0 are related to the zeros of
Px[r]: we assume that Px_r does not have a zero at x = 0. In this case series
solution to P(z,θ)u(z) = 0 can be written down as

  u(z) = sum_{0≤i<∞, 0≤j≤K} u_{i,j} z^(λ+i)*log^j(z)/j!

for some K, and which u_{i,j} are free and which u_{i,j} determined by pervious
terms in the solution is explained in Corollary 5.4: Only the u_{i,j} for
which j is strictly less than the multiplicy of λ+i as a root of Pθ_0 are free.
Considering all λ such that λ is a root of Pθ_0 and none of λ-1, λ-2, ... are
roots of Pθ[0] gives a total of deg(Pθ_0) linearly independent solutions.

Now consider a truncation

  u_N(z) = sum_{0≤i<N, 0≤j≤K} u_{i,j} z^(λ+i)*log^j(z)/j!

and the normalized difference y(z) = Px_r(z)*(u_N(z) - u(z)).
This y(z) satisfies L(z,θ)y(z) = Q_0(θ)q(z) where

  1. L(x,θ) := P(z,θ)/Px_r(z) = Q_0(θ) + Q_1(θ)*z + Q_2(θ)*z^2 + ...
  2. Q_0(θ) is monic in θ with degree r (and is proportional to Pθ_0)
  3. The Q_1, Q_2, ... all have degree < r
  4. Assuming none of λ+N, ..., λ+N+s-1 are roots of Q_0(θ), the normalized
     residue q(z) is of the form

        q(z) = sum_{0≤i<s, 0≤j≤K} q_{i,j} z^(λ+N+i)*log^j(z)/j!

Let ahat(z) be a series satisfying Proposition 5.5. Compute this by noting

    sum_{1≤j} Q_j(θ)*x^j = P(x,θ)/Px_r(x) - Q_0(θ)
                         = P(x,θ)/Px_r(x) - P(0,θ)/Px_r(0)

If we can expand the rhs as a finite linear combination

    sum_i f_i(θ)*(some power series in RR_+[[x]])

then it suffices to bound each of the f_j(θ)/Q_0(θ) in accordance with the lhs
of 5.9. Since we are ultimately interested in hhat(z) = exp(int_0^z ahat(z)/z dz),
it makes since to expand

    sum_{1≤j} Q_j(θ)*x^j = sum_i f_i(θ) x*d/dx(some power series in x*RR_+[[x]])

For all the differential equations considered here, since we have few
singularities, the "some power series in x*RR_+[x]" can be taken to be

  log(1/(1-x)),  log(1/(1-x^2)),  log((1+x)/(1-x)),
  x/(1-x)^i,     x^1/(1-x^2)^i,   x^2/(1-x^2)^i        1≤i≤~s

and our hhat(z) will take the shape (1-z)^?*(1+z)^?*exp(rational function of z).

At this point, with some minor technicalities, Algorithm 6.11 computes a ghat
such that

(*)   z^-λ*(u_N(z) - u(z)) << phat(z)*z^N*ghat(z)*hhat(z)

where phat(z) is a majorant of 1/Px_r(z). Note that we have taken out the z^N
from the ghat(z). Since the lhs is a polynomial in log(z), this must be
interpreted as saying that the rhs majorizes each of the coefficients of log(z).

Given (*) how to compute bounds on |u_N(z) - u(z)| and the derivatives
|u_N'(z) - u'(z)|, |u_N''(z) - u''(z)|, ... |u_N^(δ-1)(z) - u^(δ-1)(z)| ???

The code currently computes the power series to order O(ε^δ) of

  phat(|z|+ε)*(|z|+ε)^(λ+N)*ghat(|z|+ε)*hhat(|z|+ε)*(sum_{j<?} log^j(|z|+ε)/j!

and does an add_error! on the coefficients of ε^d/d!, but this looks suspicious.
=#

#=
Example.

Consider computing 2F1(a,b,c,z) following the "compute_f_anywhere" path.

The equation satisfied by u(x) = (1+x)^(-2*a)*2F1(a,b,c,4*x/(1+x)^2) is

  julia> R,(a,b,c)=PolynomialRing(ZZ,["a","b","c"])

  julia> F=FractionField(R)

  julia> (a,b,c)=map(F,(a,b,c))

  julia> Fx, x = PolynomialRing(F, "x")

  julia> hyp_equ([a,b], [c], -2*Fx(a)//(1+x), 4*x//(1+x)^2)

  (AbstractAlgebra.Generic.Poly{AbstractAlgebra.Generic.Frac{fmpz_mpoly}}[
    θ^2 + (c - 1)*θ,
    (-4*b + 2*c)*θ - 4*a*b + 2*a*c + 4*b - 2*c,
    -θ^2 + (-4*a + c + 3)*θ - 4*a^2 + 2*a*c + 6*a - 2*c - 2],

  AbstractAlgebra.Generic.Poly{AbstractAlgebra.Generic.Frac{fmpz_mpoly}}[
    (-4*a^2 + 2*a*c + 6*a - 2*c - 2)*x^2 + (-4*a*b + 2*a*c + 4*b - 2*c)*x,
    (-4*a + c + 3)*x^2 + (-4*b + 2*c)*x + (c - 1),
    -x^2 + 1])

which indicates

P(x, θ) = -(θ+2*a-2)*(θ+2*a-c-1)*x^2 - 2*(2*b-c)*(θ+a-1)*x + θ*(θ+c-1)
        = θ^2*(1-x)*(1+x) + ...


Now the partial fraction decomposition after integration is

  julia> Fθ,θ=PolynomialRing(F,"θ")

  julia> Fθz,z=PolynomialRing(Fθ,"z")

  julia> partial_fractions(divexact(-(θ+2*a-2)*(θ+2*a-c-1)*z^2 -
             2*(2*b-c)*(θ+a-1)*z + θ*(θ+c-1) - θ*(θ+c-1)*(1-z)*(1+z), z), 1, 1)

  ((-2*a + c + 1)*θ - 2*a^2 + a*c + 3*a - c - 1)*log(1/(1 - z^2)) +
     ((-2*b + c)*θ - 2*a*b + a*c + 2*b - c)*log((1 + z)/(1 - z))


Since Q_0(θ) = θ*(θ+c-1), For large N, our hhat(z) is going to be about

 hhat(z) = (1/(1-z^2))^|-2*a+c+1| * ((1+z)/(1-z))^|-2*b+c|,

and phat(z) is 1/(1-z^2).


From P(x, θ), the coefficients of u(z) = sum_n u_n*z^n satisfy

  (n+1)*(n+c)*u_{n+1} = 2*(2*b-c)*(n+a)*u_n + (n+2*a-1)*(n+2*a-c)*u_{n-1}.

A tight asymptotic bound on the solutions to this is

    u_n = O(n^(|2*b-c|+2*a-c-1))

so how far is the general method from asymptotically optimal?

=#

# for all differential equations here the int_0^z ahat(z)/z dz series can
# expressed as a finite linear combination of the following functions
# for vector entries, the part [i] gives the coefficient
mutable struct hyp_majorant{S}
  log_coeff::S              # log(1/(1-z))
  log2_coeff0::S            # log(1/(1-z^2))
  log2_coeff1::S            # log((1+z)/(1-z))
  poly_coeffs::Vector{S}    # [z^i]         # not used ?
  recp_coeffs::Vector{S}    # [z/(1-z)^i]
  recp2_coeffs1::Vector{S}  # [z^1/(1-z^2)^i]
  recp2_coeffs2::Vector{S}  # [z^2/(1-z^2)^i]
end

function expressify(a::hyp_majorant; context = context)
  s = Expr(:call, :+)
  push!(s.args, Expr(:call, :*, expressify(a.log_coeff), :(log(1/(1-z)))))
  push!(s.args, Expr(:call, :*, expressify(a.log2_coeff0), :(log(1/(1-z^2)))))
  push!(s.args, Expr(:call, :*, expressify(a.log2_coeff1), :(log((1+z)/(1-z)))))
  for i in 1:length(a.recp_coeffs)
    push!(s.args, Expr(:call, :*, expressify(a.recp_coeffs[i]), :(z/(1-z)^$i)))    
  end
  for i in 1:length(a.recp2_coeffs1)
    push!(s.args, Expr(:call, :*, expressify(a.recp2_coeffs1[i]), :(z/(1-z^2)^$i)))    
  end
  for i in 1:length(a.recp2_coeffs2)
    push!(s.args, Expr(:call, :*, expressify(a.recp2_coeffs1[i]), :(z^2/(1-z^2)^$i)))    
  end
  for i in 1:length(a.poly_coeffs)
    push!(s.args, Expr(:call, :*, expressify(a.poly_coeffs[i]), :(z^$i)))    
  end
  return s
end

function coeffs(a::hyp_majorant)
  return Iterators.flatten(((a.log_coeff,),
                            (a.log2_coeff0,),
                            (a.log2_coeff1,),
                            a.poly_coeffs,
                            a.recp_coeffs,
                            a.recp2_coeffs1,
                            a.recp2_coeffs2))
end

function hyp_majorant{narb}(a::hyp_majorant{S}) where S
  return hyp_majorant{narb}(
          zero(narb),
          zero(narb),
          zero(narb),
          [zero(narb) for i in 1:length(a.poly_coeffs)],
          [zero(narb) for i in 1:length(a.recp_coeffs)],
          [zero(narb) for i in 1:length(a.recp2_coeffs1)],
          [zero(narb) for i in 1:length(a.recp2_coeffs2)])
end

function Base.show(io::IO, a::hyp_majorant)
  print(io, Nemo.AbstractAlgebra.obj_to_string(a))
end

# series expansion to O(z^ord)
function series(a::hyp_majorant{narb}, ord::Int, p::Int)
  s = zero(narb_poly)
  for j in 1:ord-1
    c = zero(nacb)

#  log_coeff::S              # log(1/(1-z))
    add!(c, c, mul(a.log_coeff, fmpq(1//j), p), p)

#  log2_coeff0::S            # log(1/(1-z^2))
#  log2_coeff1::S            # log((1+z)/(1-z))
    add!(c, c, mul(iseven(j) ? a.log2_coeff0 : a.log2_coeff1, fmpq(2//j), p), p)

#  poly_coeffs::Vector{S}    # [z^i]         # not used ?
    for i in 1:length(a.poly_coeffs)
      add!(c, c, a.poly_coeffs[i], p)
    end

#  recp_coeffs::Vector{S}    # [z/(1-z)^i]
    for i in 1:length(a.recp_coeffs)
      d = binomial(ZZ(i-1+j-1), ZZ(j-1))
      add!(c, c, mul(a.poly_coeffs[i], d, p), p)      
    end

    if isodd(j)
#  recp2_coeffs1::Vector{S}  # [z^1/(1-z^2)^i]
      for i in 1:length(a.recp2_coeffs1)
        d = binomial(ZZ(i-1+div(j-1,2)), ZZ(div(j-1,2)))
        add!(c, c, mul(a.recp2_coeffs1[i], d, p), p)
      end
    else
#  recp2_coeffs2::Vector{S}  # [z^2/(1-z^2)^i]
      for i in 1:length(a.recp2_coeffs1)
        d = binomial(ZZ(i-1+div(j-2,2)), ZZ(div(j-2,2)))
        add!(c, c, mul(a.recp2_coeffs2[i], d, p), p)
      end
    end

    setcoeff!(s, j, c)
  end
  return s
end

# evaluate it at z = z
function eval_series(a::hyp_majorant{narb}, z::narb_poly, ord::Int, p::Int)
  l1 = neg!(log1p_series(-z, ord, p))
  l2 = log1p_series(z, ord, p)

#  log_coeff::S              # log(1/(1-z))
  s = mul(l1, a.log_coeff, p)

#  log2_coeff0::S            # log(1/(1-z^2))   = l1 - l2
#  log2_coeff1::S            # log((1+z)/(1-z)) = l1 + l2
  add!(s, s, mul(l1, add(a.log2_coeff1, a.log2_coeff0, p), p), p)
  add!(s, s, mul(l2, sub(a.log2_coeff1, a.log2_coeff0, p), p), p)

#  poly_coeffs::Vector{S}    # [z^i]         # not used ?
  v = z
  t = one(narb_poly)
  for i in 1:length(a.poly_coeffs)
    t = mullow(t, v, ord, p)
    add!(s, s, mul(t, a.poly_coeffs[i], p), p)
  end

#  recp_coeffs::Vector{S}    # [z/(1-z)^i]
  v = inv_series(sub(1, z, p), ord, p)
  t = z
  for i in 1:length(a.recp_coeffs)
    t = mul_low(t, t, v, ord, p)
    add!(s, s, mul(t, a.recp_coeffs[i], p), p)
  end

#  recp2_coeffs1::Vector{S}  # [z^1/(1-z^2)^i]
  v = inv_series(sub(1, mullow(z, z, ord, p), p), ord, p)
  t = z
  for i in 1:length(a.recp2_coeffs1)
    t = mullow(t, t, v, ord, p)
    add!(s, s, mul(t, a.recp2_coeffs1[i], p), p)
  end

#  recp2_coeffs2::Vector{S}  # [z^2/(1-z^2)^i]
  t = mullow(z, z, ord, p)
  for i in 1:length(a.recp2_coeffs1)
    t = mullow(t, t, v, ord, p)
    add!(s, s, mul(t, a.recp2_coeffs1[i], p), p)
  end

  return s
end

# p += a*x^i  in F[x]
function add_coeff!(p, i, a, F)
  if iszero(a)
    return
  end
  while length(p) <= i
    push!(p, zero(F))
  end
  p[1+i] += a
end

# expansion of f/((1-z)^k1*(1+z)^k2)
function partial_fractions(f, k1::Int, k2::Int)
  Fθz = parent(f)
  Fθ = base_ring(Fθz)
  F = base_ring(Fθ)
  S = elem_type(Fθ)
  z = gen(Fθz)

  log_coeff = zero(Fθ)
  log2_coeff0 = zero(Fθ)
  log2_coeff1 = zero(Fθ)
  poly_coeffs = S[]
  recp_coeffs = S[]
  recp2_coeffs1 = S[]
  recp2_coeffs2 = S[]

  if degree(f) >= k1 + k2
    q, f = divrem(f, (1-z)^k1*(1+z)^k2)
    for i in 0:degree(q)
      add_term!(poly_coeffs, i+1, divexact(coeff(q,i), 1+i), Fθ)
    end
  end

  while k1>0 || k2>0
    if k1>k2
      a = F(2)^-k2*evaluate(f, Fθ(1))
      if k1 == 1
        f = divexact(-a*(1+z)^k2+f, 1-z)
        log_coeff += a
      else
        a = divexact(a, k1-1)
        f = divexact(-a*(1+z)^k2*(1+(k1-2)*z)+f, 1-z)
        add_coeff!(recp_coeffs, k1-1, a, Fθ)
      end
      k1 -= 1
	  else
      e1 = F(2)^(k2-k1)*evaluate(f, Fθ(-1))
      f = (1-z)^(k2-k1)*f
      e2 = evaluate(f, Fθ(1))
      a = divexact(e2+e1, k2 == 1 ? 4 : 4*k2-4);
      b = divexact(e2-e1, k2 == 1 ? 4 : 4*k2-4);
      if k2 == 1
        f = divexact(-2*(a+b*z)+(1-z)^(1-k1)*f, 1-z^2)
        log2_coeff1 += a
        log2_coeff0 += b
      else
        f = divexact(a*(-1+(3-2*k2)*z^2)-2*b*(z+(-2+k2)*z^3)+f,1-z^2)
        add_coeff!(recp2_coeffs1, k2-1, a, Fθ)
        add_coeff!(recp2_coeffs2, k2-1, b, Fθ)
      end
      k1 = k2 = k2 - 1
    end
  end

  return hyp_majorant{S}(log_coeff,
                           log2_coeff0,
                           log2_coeff1,
                           poly_coeffs,
                           recp_coeffs,
                           recp2_coeffs1,
                           recp2_coeffs2)
end

# bound n*f(n)/(g(n)=prod_α(n-α)) for n0 <= n <= n1
# it is good if g does not vanish at any integer in [n0, n1]
# pass n1 < 0 for n1 = ∞ :(
function fraction_bound_normal(
  f, g,
  αs::Vector{T},
  τ::Int,
  n0::fmpz, n1::fmpz,
  p::Int) where T

  d = length(αs)
  @assert d == degree(g)

  if iszero(f)
    return zero(narb)
  end

println("fraction_bound_normal called")
@show (f, g)
@show αs
@show τ
@show (n0, n1)

  # x = interval containing [1/n1, 1/n0]
  pm_one = narb()
  ccall((:arb_zero_pm_one, libarb), Nothing,
        (Ref{narb},),
        pm_one)
  x = nacb(mul(pm_one, 1//n0, p))

  num = zero(nacb_poly)
  t = one(nacb)
  t2 = nacb_poly(one(nacb), x)
  c0 = nacb()
  c1 = nacb()
  for i in d-1:-1:0
    mul!(t, t, x, p)
    mul!(c0, t, (i>0 ? coeff(f,i-1) : 0) - coeff(f,d-1)*coeff(g,i), p)
    mul!(c1, t, -coeff(f, i), p)
    add!(num, nacb_poly(c0, c1), mullow(t2, num, τ, p), p)
  end

  den = one(nacb_poly)
  c = one(narb)
  for α in αs
    mul!(c0, x, α, p)
    sub!(c0, 1, c0, p)
    den = mullow(nacb_poly(c0, x), den, τ, p)
    if α > 0 # TODO
      loc = numerator(ceil(α))  # TODO
      if loc <= n0
        m = abs(1-α/n0)
      elseif n1 >= 0 && n1 < loc
        m = abs(1-α/n1)
      else
        m = min(abs(1-α/(loc-1)), abs(1-α/loc))
      end
      mul!(c, c, m)
    end
  end
  inv!(c, c, p)

  # make c centered around 0
  mul!(c, c, pm_one, p)

  mul!(den, den, nacb(c), p)
  setcoeff!(den, 0, one(nacb))
  rat = div_series(num, den, τ, p)
  res = zero(narb)
  for i in 0:τ-1
    add!(res, res, abs(coeff(rat, i), p), p)
  end
  mul!(res, res, c, p)
  add!(res, res, abs(nacb(coeff(f,d-1), p), p), p)

println("fraction_bound_normal returning ", res)

  return res
end

function fraction_bound_special(f, g, μ::Int, τ::Int, n0::fmpz, p::Int)
  Fθ = parent(f)
  @assert Fθ == parent(g)
  θ = gen(Fθ)
  f = evaluate(f, θ+n)
  g = evaluate(g, θ+n)
  F = nacb_poly()
  G = nacb_poly()
  for i in 0:τ-1
    setcoeff!(F, i, nacb(coeff(f, i), p))
    setcoeff!(G, i, nacb(coeff(g, μ+i), p))
  end
  rat = div(num, den, τ, p)
  res = zero(narb)
  for i in 0:τ-1
    add!(res, res, abs(coeff(rat, i)), p)
  end
  return res  
end


function equ_bound(
  Pθ, Px,
  λ::T,
  ns::Vector{fmpz},
  αs::Vector{T},
  τ::Int,
  N::Int,
  p::Int) where T

println("equ_bound called")
@show Pθ
@show Px

  Fθ = parent(Pθ[1])
  θ = gen(Fθ)
  s = length(Pθ)-1
  r = length(Px)-1
  c = coeff(Px[1+r], 0)

  # Q0 = Q_0(θ) shifted by lambda
  Q0 = evaluate(divexact(Pθ[1+0],c), θ + λ)
@show Q0

  # f = (Pθ - Px[1+r]*Q0)/(c*x) convert to bivariate with z a.k.a x on the outside
  # f is also shifted by λ
  Fθz,z = PolynomialRing(Fθ, "z")
  f = sum(evaluate(Pθ[1+i], θ + λ)*z^i for i in 0:s)
  f = divexact(f - Px[1+r](z)*Q0, c*z)
@show f

  # also shift the roots αs of Q0
  αs = αs .- λ
  @assert Q0 == prod(θ - α for α in αs)

  # denominator should be of the form (1-x)^k1*(1+x)^k2
  x = gen(parent(Px[1+r]))
  den = divexact(Px[1+r], c)
  k1 = k2 = 0
  while degree(den) > 0
    Q, R = divrem(den, 1-x)
    if iszero(R)
      k1 += 1
      den = Q
    else
      Q, R = divrem(den, 1+x)
      if iszero(R)
        k2 += 1
        den = Q
      else
        @assert false
      end
    end
  end
  @assert isone(den)
@show (k1, k2)

  pf = partial_fractions(f, k1, k2)
@show pf

  # collect the multiplicities of the ns
  nμ = Tuple{fmpz, Int}[]
  for n in ns
    if !empty(nμ) && nμ[end][1] == n
      nμ[end] = (n, nμ[end][2]+1)
    else
      push!(nμ, (n, 1))
    end
  end

  # start with 0's and take max's along intervals all the way from N to ∞
  maj = hyp_majorant{narb}(pf)
  curN = fmpz(N)
  curτ = τ
  for (n, μ) in nμ
    majorant_bound_normal(maj, pf, Q0, αs, curτ, curN, n - 1, p)
    curτ += 1
    majorant_bound_special(maj, pf, Q0, curτ, μ, curN, n, p)
    curN = fmpz(n + 1)
  end
  majorant_bound_normal(maj, pf, Q0, αs, curτ, curN, fmpz(-1), p)
@show maj
  # 1/Px[1+r] is majorized by abs(1/c)/((1-x)^max(0,k1-k2)*(1-x^2)^k2)
  return (abs(nacb(inv(c), p), p), max(0,k1-k2), k2, maj)
end

function majorant_bound_normal(
  maj::hyp_majorant{narb},
  pf::hyp_majorant{S},
  Q0::S,
  αs::T,
  τ::Int,
  n0::fmpz, n1::fmpz,
  p::Int) where {S, T}

println("majorant_bound_normal called")
@show pf
@show Q0

  for (i, j) in zip(coeffs(maj), coeffs(pf))
    max!(i, i, fraction_bound_normal(j, Q0, αs, τ, n0, n1, p), p)
  end
end

# u is a sxτxν array of the last s coefficients
# return sxτxν array of the "normalized residues"
function q_residual(Pθ, Px, u, λ::T, N::Int, τ::Int, ν::Int) where T
  Fθ = parent(Pθ[1])
  F = base_ring(Fθ)
  θ = gen(Fθ)
  c = coeff(Px[end], 0)
  s = length(Pθ)-1
  @assert length(u) == s
  @assert length(u[1]) == τ
  @assert length(u[1][1]) == ν
  b = T[coeff(Pθ[1+i](λ+N+j+θ), k) for i in 0:s, j in 0:s, k in 0:τ-1]
  q = T[F() for i in 1:s, j in 1:τ, l in 1:ν]
  v = T[F() for i in 1:s, j in 1:τ, l in 1:ν]
  for j in 0:s-1, k in τ-1:-1:0, l in 1:ν
    v[1+j,1+k,l] = zero(F)
    for jp in 1:s-j, kp in 0:τ-1-k
      v[1+j,1+k,l] += b[1+j+jp,1+j,1+kp]*u[jp][1+k+kp][l]
    end
    q[1+j,1+k,l] = c*v[1+j,1+k,l]
    for kp in 1:τ-1-k
      q[1+j,1+k,l] -= b[1+0,1+j,1+kp]*q[1+j,1+k+kp,l]
    end
    q[1+j,1+k,l] //= b[1+0,1+j,1+0]
  end
  return q
end

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
  Pθ::Vector,
  Px::Vector,
  λ::T, ns::Vector{fmpz}, # indicial roots λ+ns[1], λ+ns[2], ..., λ+ns[ν]
  αs::Vector{T},  # all indicial roots
  z::nacb,
  δ::Int,         # calculate f(z), f'(z), ... f^(δ-1)(z)
  normal::Bool,   # continuity along z < 0 is determined by counterclockwise direction
  invz::nacb,     # 1/z only used if normal = false
  p::Int) where T

println("eval_basis")
@show Pθ
@show Px
@show (λ, ns)
@show αs
@show z

  ν = length(ns)    # number of initial conditions
  τ = 0             # strict bound on the power of log(z) thus far
  maxN = 10+20*p    # max number of terms to sum

  @assert length(ns) > 0
  @assert ns[1] == 0

  s = length(Pθ) - 1
  @assert s > 0

  Fθ = parent(Pθ[1])
  θ = gen(Fθ)
  F = base_ring(Fθ)

  # u is an array of the last s solutions
  # each solution is an array of coefficients of log^k(z)/k!  0 <= k < τ
  # each coefficient is an array of ν elements of F
  u = [Vector{T}[] for i in 1:s]

  # Σ is an array of terms for the i^th derivative
  # each term is an array of coefficients of log^k(z)/k!  0 <= k < τ
  # each coefficient is an array of ν elements of nacb
  Σ = [Vector{nacb}[] for i in 1:δ]

  #
  ns = deepcopy(ns)
  iv = [[i == j ? one(F) : zero(F) for i in 1:ν] for j in 1:ν]

  t = nacb()
  zn = one(nacb)

  changedcount = 0
  N = 0
  while N < maxN
    # evaluate the next coefficient u_n
    # the evaluation is currently done in the field F and subject to blowup
    # TODO switch to ball arithmetic after a few steps
    Pn = map(a -> a(θ + (λ + N)), Pθ)
    un = get_multiplicity!(ns, iv, N)
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
    # μ = number of initial conditions at λ + N = multiplicity of root
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

    # add un*z^N to sum
    changed = false
    znrow = [[mul(zn, un[1+i][k], p) for k in 1:ν] for i in 0:τ-1]
    for d in 0:δ-1
      if d > 0
        # differentiate un
        for i in 0:τ-1, k in 1:ν          
          mul!(znrow[1+i][k], znrow[1+i][k], λ+N-d+1, p)
          if i+1 < τ
            add!(znrow[1+i][k], znrow[1+i][k], znrow[1+i+1][k], p)
          end
        end
      end
      while length(Σ[1+d]) < τ
        push!(Σ[1+d], [zero(nacb) for k in 1:ν])
      end
      for i in 0:τ-1, k in 1:ν
        add!(t, Σ[1+d][1+i][k], znrow[1+i][k], p)
        changed = changed || !overlaps(Σ[1+d][1+i][k], t)
        swap!(Σ[1+d][1+i][k], t)
      end
    end

    mul!(zn, zn, z, p)
    N += 1
    if !changed
      if (changedcount += 1) > 1
        break
      end
    else
        changedcount = 0
    end
  end

  # c/((1-z)^k1*(1-z^2)^k2) is supposed to majorize 1/Px[1+r], r = deg_θ(P)
  # exp(maj(z)) is the hhat(z)
  (c, k1, k2, maj) = equ_bound(Pθ, Px, λ, ns, αs, τ, N, p)
@show c
@show k1
@show k2
@show maj

  q = q_residual(Pθ, Px, u, λ, N, τ, ν)
@show q

  # The error bounds need more work
  finalerror = [narb_poly() for l in 1:ν]
  for l in 1:ν
    f = zero(narb_poly)
    for i in 0:s-1
      # first take the max's of coeffs of log(z)^j/j!
      m = zero(narb)
      for j in 0:s-1
        max!(m, m, abs(nacb((N+i)*q[1+i,1+j,l], p), p), p)
      end
      setcoeff!(f, i, m)
    end
    g = mullow(f, exp_series(neg!(series(maj, s, p)), s, p), s, p)
    for i in 0:s-1
      setcoeff!(g, i, max(zero(narb), div(coeff(g, i), N+i, p), p))
    end

    # f(z) is given by sum_{i,j} u_{i,j}*z^(λ+i)*log(z)^j/j!
    # The sum on j is finite and we have evaluated the sum f_N(z) on i<N.
    # For each j, the remainder sum_{i>=N} u_{i,j}*z^(λ+i) is majorized by
    #     R(z) = z^(λ+N)*g(z)*exp(maj(z))*c/((1-z)^k1*(1-z^2)^k2)
    #         TODO!!  when there are large initial roots so that ns is still
    #                 not empty increase g(z) accordingly

    # The m^th derivative error |f^m(z) - f_N^m(z)| can be bounded by the
    # ε coefficients of R(|z|+ε)*sum_{j<?}log(|z|+ε)^j/j!   ??????

    # z^(λ+N)
    zeval = narb_poly(abs(z, p), 1)
    finalerror[l] = pow_series(zeval, narb(λ + N, p), δ, p) # TODO

    # c/((1-z)^k1*(1-z^2)^k2)
    mul!(finalerror[l], finalerror[l], c, p)
    f = pow_series(inv_series(sub(1, zeval, p), δ, p), narb(k1), δ, p)
    finalerror[l] = mullow(finalerror[l], f, δ, p)
    f = pow_series(inv_series(add(1, zeval, p), δ, p), narb(k2), δ, p)
    finalerror[l] = mullow(finalerror[l], f, δ, p)

    # g(z)
    t = one(narb_poly)
    f = narb_poly(coeff(g, 0))
    for i in 1:s-1
      t = mullow(t, zeval, δ, p)
      add!(f, f, mul(t, coeff(g, i), p))
    end
    finalerror[l] = mullow(finalerror[l], f, δ, p)

    # exp(maj(z))
    f = exp_series(eval_series(maj, zeval, δ, p), δ, p)
    finalerror[l] = mullow(finalerror[l], f, δ, p)

    # sum_{j<τ} log(z)^j/j!   TODO!!! τ is not correct here if ns is not empty
    logzeval = log_series(zeval, δ, p)
    f = one(narb_poly)
    for j in τ-1:-1:1
      div!(f, mullow(f, logzeval, δ, p), narb(j, p), p)
      add!(f, f, 1, p)
    end
    finalerror[l] = mullow(finalerror[l], f, δ, p)
  end

  # evaluate the polynomials in log(z)
  logz = normal ? log(z, p) : neg!(log(invz, p))
  zλ = normal ? pow(z, λ, p) : pow(invz, -λ, p)
  M = nacb_mat(δ, ν)
  for d in 0:δ-1, l in 1:ν
    i = τ-1
    if i < 0
      t = zero(nacb)
    else
      t = Σ[1+d][1+i][l]
      while (i -= 1) >= 0
        div!(t, t, i + 1, p)
        mul!(t, t, logz, p)
        add!(t, t, Σ[1+d][1+i][l], p)
      end
    end
    mul!(t, t, d == 0 ? zλ : mul(zλ, pow(z, -d, p), p), p)
    add_error!(t, mul(coeff(finalerror[l], d), factorial(d, p), p))
    setindex!(M, 1+d, l, t)
  end
  return M
end


function eval_bases(Pθ, Px, ρ, αs, z, δ, normal, invz, p)
  M = eval_basis(Pθ, Px, ρ[1][1], ρ[1][2], αs, z, δ, normal, invz, p)
  for i in 2:length(ρ)
    M = hcat(M, eval_basis(Pθ, Px, ρ[i][1], ρ[i][2], αs, z, δ, normal, invz, p))
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
    mul!(t, t, bernoulli(2*n)//factorial(ZZ(2*n))*
               (-1)^(n+1)*(4^n-2)*x^(2*n-1), p) # TODO
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
    n = am[1+j]
    # μ = multiplicity of exponent λ + n
    μ = 1
    while j+μ < m && am[1+j+μ] == n
      μ += 1
    end
    # compute coeffs r[1+j], ..., r[1+j+μ-1]
    # of z^(λ+n), ..., z^(λ+n)*log^(μ-1)(z)/(μ-1)!
    for k in 0:m-1
      # t is the index of the term in the pFq series
      t = n - am[1+k]
      if t < 0
        continue
      end
      # compute the order required for the power series O(ε^ord)
      ord = 1
      for i in 0:m-1
        if i != k && am[1+i]-am[1+k]-t <= 0
          # a gamma factor has a simple pole
          ord += 1
        end
      end
      s = one(nacb_poly)
      # handle integer differences between numerator parameters
      for i in 0:m-1
        if i == k
          s2 = rising_factorial_series(nacb(λ+am[1+i], p), i, t, ord, p)
          s = mullow(s, s2, ord, p)
        else
          s2 = rgamma_series(nacb(λ+am[1+i], p), i, ord, p)
          s = mullow(s, s2, ord, p)
          d = am[1+i]-am[1+k]-t
          if d > 0
            s2 = gamma_series(nacb(d, p), i-k, ord, p)
            s = mullow(s, s2, ord, p)
          else
            s2 = rgamma_series(nacb(1-d, p), k-i, ord, p)
            s2 = mullow(s2, cscpi_series(ZZ(i-k), ord, p), ord, p)
            s = mullow(s, s2, ord, p)
            isodd(d) && neg!(s)
          end
        end
      end
      # no difference A - λ - am[1 + i] is an integer
      for A in as
        s2 = gamma_series(nacb(A-λ-am[1+k]-t, p), -k, ord, p)
        s = mullow(s, s2, ord, p)
      end
      # the difference B - λ - am[1 + i] could be an integer, but it's rgamma
      for B in b
        s2 = rgamma_series(nacb(B-λ-am[1+k]-t, p), -k, ord, p)
        s = mullow(s, s2, ord, p)
      end
      # now determine r[1+j], ..., r[1+j+μ-1]. They are already zeroed.
      for i in 0:min(ord, μ)-1
        f = 1//factorial(t) # TODO
        if i > 0
          f *= ZZ(k)^i//factorial(ZZ(i))
        end
        ff = mul(coeff(s, ord-1-i), f, p)
        if iseven(t)
          add!(r[1+j+i], r[1+j+i], ff, p)
        else
          sub!(r[1+j+i], r[1+j+i], ff, p)
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
  # 0 -> z
  Pθ, Px = hyp_equ(a, b, Fx(0)//Fx(1), x//Fx(1))
  αs = push!(F(1) .- b, F(0))
  e = eval_basis(Pθ, Px, F(0), [ZZ(0)], αs, z, 1, true, nacb(), p)
  return e[1, 1]  
end

function compute_f_at_∞(a::Vector{T}, b::Vector{T}, z::nacb, p::Int) where T
  F = parent(a[1])
  Fx, x = PolynomialRing(F, "x")
  # ∞ -> z
  Pθ, Px = hyp_equ(a, b, Fx(0)//Fx(1), -Fx(1)//x)
  αs = a
  nz = neg!(nacb(), z)
  rz = inv(nz, p)
  s = zero(nacb)
  for (λ, ns, indices) in partition_mod1(a)
    arest = T[a[j] for j in 1:length(a) if !(j in indices)]
    iv = hyp_initial_values_at_∞(λ, ns, arest, b, p)
    e = eval_basis(Pθ, Px, λ, ns, αs, rz, 1, false, nz, p)
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
  Pθ, Px = hyp_equ(a, b, Fx(0)//Fx(1), x//Fx(1))
  αs = push!(F(1) .- b, F(0))
  e0 = eval_basis(Pθ, Px, F(0), [ZZ(0)], αs, w, n, true, nacb(), p)
  for i in 2:2:n
    setindex!(e0, i, 1, neg!(e0[i, 1])) # correct odd derivatives
  end
  # roots of indicial equation are 0:n-2 and σ
  σ::T = sum(b) - sum(a)
  αs = push!(FF.(0:n-2), σ)
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
  Pθ, Px = hyp_equ(a, b, Fx(0)//Fx(1), (1-x)//Fx(1))
  e1 = eval_bases(Pθ, Px, ρ, αs, 1-z, 1, true, nacb(), p)
  e2 = eval_bases(Pθ, Px, ρ, αs, 1-w, n, true, nacb(), p)
  # 0 -> w -> 1 -> z
  return mul(e1, solve(e2, e0, p), p)[1, 1]
end

function compute_f_anywhere(a::Vector{T}, b::Vector{T}, z::nacb, p::Int) where T
  F = parent(a[1])
  Fx, x = PolynomialRing(F, "x")
  # 0 -> (1-sqrt(1-z))/(1+sqrt(1-z))
  Pθ, Px = hyp_equ(a, b, -2*Fx(a[1])//(1+x), 4*x//(1+x)^2)
  αs = push!(F(1) .- b, F(0))
  t = sqrt(1-z, p)
  s = 1+t
  e = eval_basis(Pθ, Px, F(0), [ZZ(0)], αs, div(1-t, s, p), 1, true, nacb(), p)
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
                  nacb(QQ(1//2), 100),
                  100)

#=
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
=#

nothing
