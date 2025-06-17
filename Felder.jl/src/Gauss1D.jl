#=
1D Gauss-Legendre quadrature weights and points, calculated with

https://github.com/JuliaApproximation/FastGaussQuadrature.jl
=#


##############################
# Line
##############################
#=
      ○────────────○
    -0.5          0.5

Changed parametrization from (-1, 1) to (-0.5, 0.5) to be consistent with other elements.
=#

# --------------------------- Strength φ = 1 ---------------------------
# Number of points: 1

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{1}) = 1

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{1}, i::Integer)::Float64
    i ==  1 && return 2.0000000000000000e+00 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{1} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{1}, i::Integer)::SVector{1, Float64}
    i ==  1 && return ( 0.0000000000000000e+00,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{1} in AbstractEdge"))
end

# --------------------------- Strength φ = 2 ---------------------------
# Number of points: 2

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{2}) = 2

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{2}, i::Integer)::Float64
    i == 1  && return 1.0000000000000000e+00 * 0.5
    i == 2  && return 1.0000000000000000e+00 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{2} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{2}, i::Integer)::SVector{1, Float64}
    i == 1  && return (-5.7735026918962584e-01,) .* 0.5
    i == 2  && return ( 5.7735026918962584e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{2} in AbstractEdge"))
end

# --------------------------- Strength φ = 3 ---------------------------
# Number of points: 3

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{3}) = 3

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{3}, i::Integer)::Float64
    i ==  1 && return 5.5555555555555558e-01 * 0.5
    i ==  2 && return 8.8888888888888884e-01 * 0.5
    i ==  3 && return 5.5555555555555558e-01 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{3} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{3}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-7.7459666924148340e-01,) .* 0.5
    i ==  2 && return ( 0.0000000000000000e+00,) .* 0.5
    i ==  3 && return ( 7.7459666924148340e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{3} in AbstractEdge"))
end

# --------------------------- Strength φ = 4 ---------------------------
# Number of points: 4

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{4}) = 4

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{4}, i::Integer)::Float64
    i ==  1 && return 3.4785484513745385e-01 * 0.5
    i ==  2 && return 6.5214515486254621e-01 * 0.5
    i ==  3 && return 6.5214515486254621e-01 * 0.5
    i ==  4 && return 3.4785484513745385e-01 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{4} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{4}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-8.6113631159405257e-01,) .* 0.5
    i ==  2 && return (-3.3998104358485631e-01,) .* 0.5
    i ==  3 && return ( 3.3998104358485631e-01,) .* 0.5
    i ==  4 && return ( 8.6113631159405257e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{4} in AbstractEdge"))
end

# --------------------------- Strength φ = 5 ---------------------------
# Number of points: 5

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{5}) = 5

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{5}, i::Integer)::Float64
    i ==  1 && return 2.3692688505618908e-01 * 0.5
    i ==  2 && return 4.7862867049936647e-01 * 0.5
    i ==  3 && return 5.6888888888888889e-01 * 0.5
    i ==  4 && return 4.7862867049936647e-01 * 0.5
    i ==  5 && return 2.3692688505618908e-01 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{5} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{5}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.0617984593866396e-01,) .* 0.5
    i ==  2 && return (-5.3846931010568311e-01,) .* 0.5
    i ==  3 && return ( 0.0000000000000000e+00,) .* 0.5
    i ==  4 && return ( 5.3846931010568311e-01,) .* 0.5
    i ==  5 && return ( 9.0617984593866396e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{5} in AbstractEdge"))
end

# --------------------------- Strength φ = 6 ---------------------------
# Number of points: 6

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{6}) = 6

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{6}, i::Integer)::Float64
    i ==  1 && return 1.7132449237917025e-01 * 0.5
    i ==  2 && return 3.6076157304813850e-01 * 0.5
    i ==  3 && return 4.6791393457269126e-01 * 0.5
    i ==  4 && return 4.6791393457269126e-01 * 0.5
    i ==  5 && return 3.6076157304813850e-01 * 0.5
    i ==  6 && return 1.7132449237917025e-01 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{6} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{6}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.3246951420315205e-01,) .* 0.5
    i ==  2 && return (-6.6120938646626448e-01,) .* 0.5
    i ==  3 && return (-2.3861918608319690e-01,) .* 0.5
    i ==  4 && return ( 2.3861918608319690e-01,) .* 0.5
    i ==  5 && return ( 6.6120938646626448e-01,) .* 0.5
    i ==  6 && return ( 9.3246951420315205e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{6} in AbstractEdge"))
end

# --------------------------- Strength φ = 7 ---------------------------
# Number of points: 7

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{7}) = 7

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{7}, i::Integer)::Float64
    i ==  1 && return 1.2948496616887020e-01 * 0.5
    i ==  2 && return 2.7970539148927659e-01 * 0.5
    i ==  3 && return 3.8183005050511892e-01 * 0.5
    i ==  4 && return 4.1795918367346940e-01 * 0.5
    i ==  5 && return 3.8183005050511892e-01 * 0.5
    i ==  6 && return 2.7970539148927659e-01 * 0.5
    i ==  7 && return 1.2948496616887020e-01 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{7} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{7}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.4910791234275860e-01,) .* 0.5
    i ==  2 && return (-7.4153118559939446e-01,) .* 0.5
    i ==  3 && return (-4.0584515137739718e-01,) .* 0.5
    i ==  4 && return ( 0.0000000000000000e+00,) .* 0.5
    i ==  5 && return ( 4.0584515137739718e-01,) .* 0.5
    i ==  6 && return ( 7.4153118559939446e-01,) .* 0.5
    i ==  7 && return ( 9.4910791234275860e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{7} in AbstractEdge"))
end

# --------------------------- Strength φ = 8 ---------------------------
# Number of points: 8

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{8}) = 8

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{8}, i::Integer)::Float64
    i ==  1 && return 1.0122853629037676e-01 * 0.5
    i ==  2 && return 2.2238103445337445e-01 * 0.5
    i ==  3 && return 3.1370664587788744e-01 * 0.5
    i ==  4 && return 3.6268378337836193e-01 * 0.5
    i ==  5 && return 3.6268378337836193e-01 * 0.5
    i ==  6 && return 3.1370664587788744e-01 * 0.5
    i ==  7 && return 2.2238103445337445e-01 * 0.5
    i ==  8 && return 1.0122853629037676e-01 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{8} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{8}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.6028985649753629e-01,) .* 0.5
    i ==  2 && return (-7.9666647741362673e-01,) .* 0.5
    i ==  3 && return (-5.2553240991632899e-01,) .* 0.5
    i ==  4 && return (-1.8343464249564981e-01,) .* 0.5
    i ==  5 && return ( 1.8343464249564981e-01,) .* 0.5
    i ==  6 && return ( 5.2553240991632899e-01,) .* 0.5
    i ==  7 && return ( 7.9666647741362673e-01,) .* 0.5
    i ==  8 && return ( 9.6028985649753629e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{8} in AbstractEdge"))
end

# --------------------------- Strength φ = 9 ---------------------------
# Number of points: 9

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{9}) = 9

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{9}, i::Integer)::Float64
    i ==  1 && return 8.1274388361574371e-02 * 0.5
    i ==  2 && return 1.8064816069485742e-01 * 0.5
    i ==  3 && return 2.6061069640293538e-01 * 0.5
    i ==  4 && return 3.1234707704000275e-01 * 0.5
    i ==  5 && return 3.3023935500125978e-01 * 0.5
    i ==  6 && return 3.1234707704000275e-01 * 0.5
    i ==  7 && return 2.6061069640293538e-01 * 0.5
    i ==  8 && return 1.8064816069485742e-01 * 0.5
    i ==  9 && return 8.1274388361574371e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{9} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{9}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.6816023950762609e-01,) .* 0.5
    i ==  2 && return (-8.3603110732663577e-01,) .* 0.5
    i ==  3 && return (-6.1337143270059036e-01,) .* 0.5
    i ==  4 && return (-3.2425342340380892e-01,) .* 0.5
    i ==  5 && return ( 0.0000000000000000e+00,) .* 0.5
    i ==  6 && return ( 3.2425342340380892e-01,) .* 0.5
    i ==  7 && return ( 6.1337143270059036e-01,) .* 0.5
    i ==  8 && return ( 8.3603110732663577e-01,) .* 0.5
    i ==  9 && return ( 9.6816023950762609e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{9} in AbstractEdge"))
end

# --------------------------- Strength φ = 10 ---------------------------
# Number of points: 10

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{10}) = 10

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{10}, i::Integer)::Float64
    i ==  1 && return 6.6671344308688207e-02 * 0.5
    i ==  2 && return 1.4945134915058056e-01 * 0.5
    i ==  3 && return 2.1908636251598207e-01 * 0.5
    i ==  4 && return 2.6926671930999652e-01 * 0.5
    i ==  5 && return 2.9552422471475293e-01 * 0.5
    i ==  6 && return 2.9552422471475293e-01 * 0.5
    i ==  7 && return 2.6926671930999652e-01 * 0.5
    i ==  8 && return 2.1908636251598207e-01 * 0.5
    i ==  9 && return 1.4945134915058056e-01 * 0.5
    i == 10 && return 6.6671344308688207e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{10} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{10}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.7390652851717174e-01,) .* 0.5
    i ==  2 && return (-8.6506336668898454e-01,) .* 0.5
    i ==  3 && return (-6.7940956829902444e-01,) .* 0.5
    i ==  4 && return (-4.3339539412924721e-01,) .* 0.5
    i ==  5 && return (-1.4887433898163122e-01,) .* 0.5
    i ==  6 && return ( 1.4887433898163122e-01,) .* 0.5
    i ==  7 && return ( 4.3339539412924721e-01,) .* 0.5
    i ==  8 && return ( 6.7940956829902444e-01,) .* 0.5
    i ==  9 && return ( 8.6506336668898454e-01,) .* 0.5
    i == 10 && return ( 9.7390652851717174e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{10} in AbstractEdge"))
end

# --------------------------- Strength φ = 11 ---------------------------
# Number of points: 11

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{11}) = 11

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{11}, i::Integer)::Float64
    i ==  1 && return 5.5668567116173705e-02 * 0.5
    i ==  2 && return 1.2558036946490453e-01 * 0.5
    i ==  3 && return 1.8629021092773429e-01 * 0.5
    i ==  4 && return 2.3319376459199057e-01 * 0.5
    i ==  5 && return 2.6280454451024665e-01 * 0.5
    i ==  6 && return 2.7292508677790062e-01 * 0.5
    i ==  7 && return 2.6280454451024665e-01 * 0.5
    i ==  8 && return 2.3319376459199057e-01 * 0.5
    i ==  9 && return 1.8629021092773429e-01 * 0.5
    i == 10 && return 1.2558036946490453e-01 * 0.5
    i == 11 && return 5.5668567116173705e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{11} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{11}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.7822865814605697e-01,) .* 0.5
    i ==  2 && return (-8.8706259976809532e-01,) .* 0.5
    i ==  3 && return (-7.3015200557404936e-01,) .* 0.5
    i ==  4 && return (-5.1909612920681181e-01,) .* 0.5
    i ==  5 && return (-2.6954315595234496e-01,) .* 0.5
    i ==  6 && return ( 0.0000000000000000e+00,) .* 0.5
    i ==  7 && return ( 2.6954315595234496e-01,) .* 0.5
    i ==  8 && return ( 5.1909612920681181e-01,) .* 0.5
    i ==  9 && return ( 7.3015200557404936e-01,) .* 0.5
    i == 10 && return ( 8.8706259976809532e-01,) .* 0.5
    i == 11 && return ( 9.7822865814605697e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{11} in AbstractEdge"))
end

# --------------------------- Strength φ = 12 ---------------------------
# Number of points: 12

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{12}) = 12

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{12}, i::Integer)::Float64
    i ==  1 && return 4.7175336386511751e-02 * 0.5
    i ==  2 && return 1.0693932599531826e-01 * 0.5
    i ==  3 && return 1.6007832854334625e-01 * 0.5
    i ==  4 && return 2.0316742672306587e-01 * 0.5
    i ==  5 && return 2.3349253653835492e-01 * 0.5
    i ==  6 && return 2.4914704581340288e-01 * 0.5
    i ==  7 && return 2.4914704581340288e-01 * 0.5
    i ==  8 && return 2.3349253653835492e-01 * 0.5
    i ==  9 && return 2.0316742672306587e-01 * 0.5
    i == 10 && return 1.6007832854334625e-01 * 0.5
    i == 11 && return 1.0693932599531826e-01 * 0.5
    i == 12 && return 4.7175336386511751e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{12} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{12}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.8156063424671924e-01,) .* 0.5
    i ==  2 && return (-9.0411725637047480e-01,) .* 0.5
    i ==  3 && return (-7.6990267419430469e-01,) .* 0.5
    i ==  4 && return (-5.8731795428661748e-01,) .* 0.5
    i ==  5 && return (-3.6783149899818018e-01,) .* 0.5
    i ==  6 && return (-1.2523340851146891e-01,) .* 0.5
    i ==  7 && return ( 1.2523340851146891e-01,) .* 0.5
    i ==  8 && return ( 3.6783149899818018e-01,) .* 0.5
    i ==  9 && return ( 5.8731795428661748e-01,) .* 0.5
    i == 10 && return ( 7.6990267419430469e-01,) .* 0.5
    i == 11 && return ( 9.0411725637047480e-01,) .* 0.5
    i == 12 && return ( 9.8156063424671924e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{12} in AbstractEdge"))
end

# --------------------------- Strength φ = 13 ---------------------------
# Number of points: 13

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{13}) = 13

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{13}, i::Integer)::Float64
    i ==  1 && return 4.0484004765316217e-02 * 0.5
    i ==  2 && return 9.2121499837728382e-02 * 0.5
    i ==  3 && return 1.3887351021978733e-01 * 0.5
    i ==  4 && return 1.7814598076194565e-01 * 0.5
    i ==  5 && return 2.0781604753688848e-01 * 0.5
    i ==  6 && return 2.2628318026289729e-01 * 0.5
    i ==  7 && return 2.3255155323087390e-01 * 0.5
    i ==  8 && return 2.2628318026289729e-01 * 0.5
    i ==  9 && return 2.0781604753688848e-01 * 0.5
    i == 10 && return 1.7814598076194565e-01 * 0.5
    i == 11 && return 1.3887351021978733e-01 * 0.5
    i == 12 && return 9.2121499837728382e-02 * 0.5
    i == 13 && return 4.0484004765316217e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{13} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{13}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.8418305471858814e-01,) .* 0.5
    i ==  2 && return (-9.1759839922297792e-01,) .* 0.5
    i ==  3 && return (-8.0157809073330988e-01,) .* 0.5
    i ==  4 && return (-6.4234933944034023e-01,) .* 0.5
    i ==  5 && return (-4.4849275103644687e-01,) .* 0.5
    i ==  6 && return (-2.3045831595513480e-01,) .* 0.5
    i ==  7 && return ( 0.0000000000000000e+00,) .* 0.5
    i ==  8 && return ( 2.3045831595513480e-01,) .* 0.5
    i ==  9 && return ( 4.4849275103644687e-01,) .* 0.5
    i == 10 && return ( 6.4234933944034023e-01,) .* 0.5
    i == 11 && return ( 8.0157809073330988e-01,) .* 0.5
    i == 12 && return ( 9.1759839922297792e-01,) .* 0.5
    i == 13 && return ( 9.8418305471858814e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{13} in AbstractEdge"))
end

# --------------------------- Strength φ = 14 ---------------------------
# Number of points: 14

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{14}) = 14

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{14}, i::Integer)::Float64
    i ==  1 && return 3.5119460331751874e-02 * 0.5
    i ==  2 && return 8.0158087159760319e-02 * 0.5
    i ==  3 && return 1.2151857068790312e-01 * 0.5
    i ==  4 && return 1.5720316715819357e-01 * 0.5
    i ==  5 && return 1.8553839747793779e-01 * 0.5
    i ==  6 && return 2.0519846372129560e-01 * 0.5
    i ==  7 && return 2.1526385346315777e-01 * 0.5
    i ==  8 && return 2.1526385346315777e-01 * 0.5
    i ==  9 && return 2.0519846372129560e-01 * 0.5
    i == 10 && return 1.8553839747793779e-01 * 0.5
    i == 11 && return 1.5720316715819357e-01 * 0.5
    i == 12 && return 1.2151857068790312e-01 * 0.5
    i == 13 && return 8.0158087159760319e-02 * 0.5
    i == 14 && return 3.5119460331751874e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature14} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{14}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.8628380869681231e-01,) .* 0.5
    i ==  2 && return (-9.2843488366357352e-01,) .* 0.5
    i ==  3 && return (-8.2720131506976502e-01,) .* 0.5
    i ==  4 && return (-6.8729290481168548e-01,) .* 0.5
    i ==  5 && return (-5.1524863635815410e-01,) .* 0.5
    i ==  6 && return (-3.1911236892788974e-01,) .* 0.5
    i ==  7 && return (-1.0805494870734364e-01,) .* 0.5
    i ==  8 && return ( 1.0805494870734364e-01,) .* 0.5
    i ==  9 && return ( 3.1911236892788974e-01,) .* 0.5
    i == 10 && return ( 5.1524863635815410e-01,) .* 0.5
    i == 11 && return ( 6.8729290481168548e-01,) .* 0.5
    i == 12 && return ( 8.2720131506976502e-01,) .* 0.5
    i == 13 && return ( 9.2843488366357352e-01,) .* 0.5
    i == 14 && return ( 9.8628380869681231e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{14} in AbstractEdge"))
end

# --------------------------- Strength φ = 15 ---------------------------
# Number of points: 15

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{15}) = 15

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{15}, i::Integer)::Float64
    i ==  1 && return 3.0753241996117630e-02 * 0.5
    i ==  2 && return 7.0366047488108374e-02 * 0.5
    i ==  3 && return 1.0715922046717200e-01 * 0.5
    i ==  4 && return 1.3957067792615432e-01 * 0.5
    i ==  5 && return 1.6626920581699403e-01 * 0.5
    i ==  6 && return 1.8616100001556221e-01 * 0.5
    i ==  7 && return 1.9843148532711158e-01 * 0.5
    i ==  8 && return 2.0257824192556129e-01 * 0.5
    i ==  9 && return 1.9843148532711158e-01 * 0.5
    i == 10 && return 1.8616100001556221e-01 * 0.5
    i == 11 && return 1.6626920581699403e-01 * 0.5
    i == 12 && return 1.3957067792615432e-01 * 0.5
    i == 13 && return 1.0715922046717200e-01 * 0.5
    i == 14 && return 7.0366047488108374e-02 * 0.5
    i == 15 && return 3.0753241996117630e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{15} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{15}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.8799251802048549e-01,) .* 0.5
    i ==  2 && return (-9.3727339240070595e-01,) .* 0.5
    i ==  3 && return (-8.4820658341042721e-01,) .* 0.5
    i ==  4 && return (-7.2441773136017007e-01,) .* 0.5
    i ==  5 && return (-5.7097217260853883e-01,) .* 0.5
    i ==  6 && return (-3.9415134707756339e-01,) .* 0.5
    i ==  7 && return (-2.0119409399743451e-01,) .* 0.5
    i ==  8 && return ( 0.0000000000000000e+00,) .* 0.5
    i ==  9 && return ( 2.0119409399743451e-01,) .* 0.5
    i == 10 && return ( 3.9415134707756339e-01,) .* 0.5
    i == 11 && return ( 5.7097217260853883e-01,) .* 0.5
    i == 12 && return ( 7.2441773136017007e-01,) .* 0.5
    i == 13 && return ( 8.4820658341042721e-01,) .* 0.5
    i == 14 && return ( 9.3727339240070595e-01,) .* 0.5
    i == 15 && return ( 9.8799251802048549e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{15} in AbstractEdge"))
end

# --------------------------- Strength φ = 16 ---------------------------
# Number of points: 16

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{16}) = 16

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{16}, i::Integer)::Float64
    i ==  1 && return 2.7152459411754121e-02 * 0.5
    i ==  2 && return 6.2253523938647803e-02 * 0.5
    i ==  3 && return 9.5158511682492786e-02 * 0.5
    i ==  4 && return 1.2462897125553386e-01 * 0.5
    i ==  5 && return 1.4959598881657671e-01 * 0.5
    i ==  6 && return 1.6915651939500262e-01 * 0.5
    i ==  7 && return 1.8260341504492350e-01 * 0.5
    i ==  8 && return 1.8945061045506850e-01 * 0.5
    i ==  9 && return 1.8945061045506850e-01 * 0.5
    i == 10 && return 1.8260341504492350e-01 * 0.5
    i == 11 && return 1.6915651939500262e-01 * 0.5
    i == 12 && return 1.4959598881657671e-01 * 0.5
    i == 13 && return 1.2462897125553386e-01 * 0.5
    i == 14 && return 9.5158511682492786e-02 * 0.5
    i == 15 && return 6.2253523938647803e-02 * 0.5
    i == 16 && return 2.7152459411754121e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{16} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{16}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.8940093499164994e-01,) .* 0.5
    i ==  2 && return (-9.4457502307323260e-01,) .* 0.5
    i ==  3 && return (-8.6563120238783176e-01,) .* 0.5
    i ==  4 && return (-7.5540440835500300e-01,) .* 0.5
    i ==  5 && return (-6.1787624440264377e-01,) .* 0.5
    i ==  6 && return (-4.5801677765722737e-01,) .* 0.5
    i ==  7 && return (-2.8160355077925892e-01,) .* 0.5
    i ==  8 && return (-9.5012509837637441e-02,) .* 0.5
    i ==  9 && return ( 9.5012509837637441e-02,) .* 0.5
    i == 10 && return ( 2.8160355077925892e-01,) .* 0.5
    i == 11 && return ( 4.5801677765722737e-01,) .* 0.5
    i == 12 && return ( 6.1787624440264377e-01,) .* 0.5
    i == 13 && return ( 7.5540440835500300e-01,) .* 0.5
    i == 14 && return ( 8.6563120238783176e-01,) .* 0.5
    i == 15 && return ( 9.4457502307323260e-01,) .* 0.5
    i == 16 && return ( 9.8940093499164994e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{16} in AbstractEdge"))
end

# --------------------------- Strength φ = 17 ---------------------------
# Number of points: 17

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{17}) = 17

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{17}, i::Integer)::Float64
    i ==  1 && return 2.4148302868548521e-02 * 0.5
    i ==  2 && return 5.5459529373987120e-02 * 0.5
    i ==  3 && return 8.5036148317179192e-02 * 0.5
    i ==  4 && return 1.1188384719340401e-01 * 0.5
    i ==  5 && return 1.3513636846852553e-01 * 0.5
    i ==  6 && return 1.5404576107681037e-01 * 0.5
    i ==  7 && return 1.6800410215645001e-01 * 0.5
    i ==  8 && return 1.7656270536699262e-01 * 0.5
    i ==  9 && return 1.7944647035620653e-01 * 0.5
    i == 10 && return 1.7656270536699262e-01 * 0.5
    i == 11 && return 1.6800410215645001e-01 * 0.5
    i == 12 && return 1.5404576107681037e-01 * 0.5
    i == 13 && return 1.3513636846852553e-01 * 0.5
    i == 14 && return 1.1188384719340401e-01 * 0.5
    i == 15 && return 8.5036148317179192e-02 * 0.5
    i == 16 && return 5.5459529373987120e-02 * 0.5
    i == 17 && return 2.4148302868548521e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{17} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{17}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9057547531441736e-01,) .* 0.5
    i ==  2 && return (-9.5067552176876780e-01,) .* 0.5
    i ==  3 && return (-8.8023915372698591e-01,) .* 0.5
    i ==  4 && return (-7.8151400389680137e-01,) .* 0.5
    i ==  5 && return (-6.5767115921669073e-01,) .* 0.5
    i ==  6 && return (-5.1269053708647694e-01,) .* 0.5
    i ==  7 && return (-3.5123176345387630e-01,) .* 0.5
    i ==  8 && return (-1.7848418149584785e-01,) .* 0.5
    i ==  9 && return ( 0.0000000000000000e+00,) .* 0.5
    i == 10 && return ( 1.7848418149584785e-01,) .* 0.5
    i == 11 && return ( 3.5123176345387630e-01,) .* 0.5
    i == 12 && return ( 5.1269053708647694e-01,) .* 0.5
    i == 13 && return ( 6.5767115921669073e-01,) .* 0.5
    i == 14 && return ( 7.8151400389680137e-01,) .* 0.5
    i == 15 && return ( 8.8023915372698591e-01,) .* 0.5
    i == 16 && return ( 9.5067552176876780e-01,) .* 0.5
    i == 17 && return ( 9.9057547531441736e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{17} in AbstractEdge"))
end

# --------------------------- Strength φ = 18 ---------------------------
# Number of points: 18

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{18}) = 18

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{18}, i::Integer)::Float64
    i ==  1 && return 2.1616013526483541e-02 * 0.5
    i ==  2 && return 4.9714548894969846e-02 * 0.5
    i ==  3 && return 7.6425730254889052e-02 * 0.5
    i ==  4 && return 1.0094204410628715e-01 * 0.5
    i ==  5 && return 1.2255520671147846e-01 * 0.5
    i ==  6 && return 1.4064291467065065e-01 * 0.5
    i ==  7 && return 1.5468467512626519e-01 * 0.5
    i ==  8 && return 1.6427648374583262e-01 * 0.5
    i ==  9 && return 1.6914238296314354e-01 * 0.5
    i == 10 && return 1.6914238296314354e-01 * 0.5
    i == 11 && return 1.6427648374583262e-01 * 0.5
    i == 12 && return 1.5468467512626519e-01 * 0.5
    i == 13 && return 1.4064291467065065e-01 * 0.5
    i == 14 && return 1.2255520671147846e-01 * 0.5
    i == 15 && return 1.0094204410628715e-01 * 0.5
    i == 16 && return 7.6425730254889052e-02 * 0.5
    i == 17 && return 4.9714548894969846e-02 * 0.5
    i == 18 && return 2.1616013526483541e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{18} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{18}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9156516842093090e-01,) .* 0.5
    i ==  2 && return (-9.5582394957139771e-01,) .* 0.5
    i ==  3 && return (-8.9260246649755581e-01,) .* 0.5
    i ==  4 && return (-8.0370495897252314e-01,) .* 0.5
    i ==  5 && return (-6.9168704306035322e-01,) .* 0.5
    i ==  6 && return (-5.5977083107394754e-01,) .* 0.5
    i ==  7 && return (-4.1175116146284263e-01,) .* 0.5
    i ==  8 && return (-2.5188622569150554e-01,) .* 0.5
    i ==  9 && return (-8.4775013041735306e-02,) .* 0.5
    i == 10 && return ( 8.4775013041735306e-02,) .* 0.5
    i == 11 && return ( 2.5188622569150554e-01,) .* 0.5
    i == 12 && return ( 4.1175116146284263e-01,) .* 0.5
    i == 13 && return ( 5.5977083107394754e-01,) .* 0.5
    i == 14 && return ( 6.9168704306035322e-01,) .* 0.5
    i == 15 && return ( 8.0370495897252314e-01,) .* 0.5
    i == 16 && return ( 8.9260246649755581e-01,) .* 0.5
    i == 17 && return ( 9.5582394957139771e-01,) .* 0.5
    i == 18 && return ( 9.9156516842093090e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{18} in AbstractEdge"))
end

# --------------------------- Strength φ = 19 ---------------------------
# Number of points: 19

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{19}) = 19

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{19}, i::Integer)::Float64
    i ==  1 && return 1.9461788229725972e-02 * 0.5
    i ==  2 && return 4.4814226765699586e-02 * 0.5
    i ==  3 && return 6.9044542737641115e-02 * 0.5
    i ==  4 && return 9.1490021622449985e-02 * 0.5
    i ==  5 && return 1.1156664554733405e-01 * 0.5
    i ==  6 && return 1.2875396253933630e-01 * 0.5
    i ==  7 && return 1.4260670217360646e-01 * 0.5
    i ==  8 && return 1.5276604206585964e-01 * 0.5
    i ==  9 && return 1.5896884339395442e-01 * 0.5
    i == 10 && return 1.6105444984878370e-01 * 0.5
    i == 11 && return 1.5896884339395442e-01 * 0.5
    i == 12 && return 1.5276604206585964e-01 * 0.5
    i == 13 && return 1.4260670217360646e-01 * 0.5
    i == 14 && return 1.2875396253933630e-01 * 0.5
    i == 15 && return 1.1156664554733405e-01 * 0.5
    i == 16 && return 9.1490021622449985e-02 * 0.5
    i == 17 && return 6.9044542737641115e-02 * 0.5
    i == 18 && return 4.4814226765699586e-02 * 0.5
    i == 19 && return 1.9461788229725972e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{19} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{19}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9240684384358435e-01,) .* 0.5
    i ==  2 && return (-9.6020815213483002e-01,) .* 0.5
    i ==  3 && return (-9.0315590361481790e-01,) .* 0.5
    i ==  4 && return (-8.2271465653714282e-01,) .* 0.5
    i ==  5 && return (-7.2096617733522939e-01,) .* 0.5
    i ==  6 && return (-6.0054530466168099e-01,) .* 0.5
    i ==  7 && return (-4.6457074137596094e-01,) .* 0.5
    i ==  8 && return (-3.1656409996362983e-01,) .* 0.5
    i ==  9 && return (-1.6035864564022537e-01,) .* 0.5
    i == 10 && return ( 0.0000000000000000e+00,) .* 0.5
    i == 11 && return ( 1.6035864564022537e-01,) .* 0.5
    i == 12 && return ( 3.1656409996362983e-01,) .* 0.5
    i == 13 && return ( 4.6457074137596094e-01,) .* 0.5
    i == 14 && return ( 6.0054530466168099e-01,) .* 0.5
    i == 15 && return ( 7.2096617733522939e-01,) .* 0.5
    i == 16 && return ( 8.2271465653714282e-01,) .* 0.5
    i == 17 && return ( 9.0315590361481790e-01,) .* 0.5
    i == 18 && return ( 9.6020815213483002e-01,) .* 0.5
    i == 19 && return ( 9.9240684384358435e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{19} in AbstractEdge"))
end

# --------------------------- Strength φ = 20 ---------------------------
# Number of points: 20

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{20}) = 20

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{20}, i::Integer)::Float64
    i ==  1 && return 1.7614007139152055e-02 * 0.5
    i ==  2 && return 4.0601429800386973e-02 * 0.5
    i ==  3 && return 6.2672048334109040e-02 * 0.5
    i ==  4 && return 8.3276741576704574e-02 * 0.5
    i ==  5 && return 1.0193011981724034e-01 * 0.5
    i ==  6 && return 1.1819453196151841e-01 * 0.5
    i ==  7 && return 1.3168863844917655e-01 * 0.5
    i ==  8 && return 1.4209610931838201e-01 * 0.5
    i ==  9 && return 1.4917298647260377e-01 * 0.5
    i == 10 && return 1.5275338713072587e-01 * 0.5
    i == 11 && return 1.5275338713072587e-01 * 0.5
    i == 12 && return 1.4917298647260377e-01 * 0.5
    i == 13 && return 1.4209610931838201e-01 * 0.5
    i == 14 && return 1.3168863844917655e-01 * 0.5
    i == 15 && return 1.1819453196151841e-01 * 0.5
    i == 16 && return 1.0193011981724034e-01 * 0.5
    i == 17 && return 8.3276741576704574e-02 * 0.5
    i == 18 && return 6.2672048334109040e-02 * 0.5
    i == 19 && return 4.0601429800386973e-02 * 0.5
    i == 20 && return 1.7614007139152055e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{20} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{20}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9312859918509488e-01,) .* 0.5
    i ==  2 && return (-9.6397192727791381e-01,) .* 0.5
    i ==  3 && return (-9.1223442825132595e-01,) .* 0.5
    i ==  4 && return (-8.3911697182221878e-01,) .* 0.5
    i ==  5 && return (-7.4633190646015080e-01,) .* 0.5
    i ==  6 && return (-6.3605368072651502e-01,) .* 0.5
    i ==  7 && return (-5.1086700195082713e-01,) .* 0.5
    i ==  8 && return (-3.7370608871541955e-01,) .* 0.5
    i ==  9 && return (-2.2778585114164510e-01,) .* 0.5
    i == 10 && return (-7.6526521133497324e-02,) .* 0.5
    i == 11 && return ( 7.6526521133497324e-02,) .* 0.5
    i == 12 && return ( 2.2778585114164510e-01,) .* 0.5
    i == 13 && return ( 3.7370608871541955e-01,) .* 0.5
    i == 14 && return ( 5.1086700195082713e-01,) .* 0.5
    i == 15 && return ( 6.3605368072651502e-01,) .* 0.5
    i == 16 && return ( 7.4633190646015080e-01,) .* 0.5
    i == 17 && return ( 8.3911697182221878e-01,) .* 0.5
    i == 18 && return ( 9.1223442825132595e-01,) .* 0.5
    i == 19 && return ( 9.6397192727791381e-01,) .* 0.5
    i == 20 && return ( 9.9312859918509488e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{20} in AbstractEdge"))
end

# --------------------------- Strength φ = 21 ---------------------------
# Number of points: 21

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{21}) = 21

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{21}, i::Integer)::Float64
    i ==  1 && return 1.6017228257773849e-02 * 0.5
    i ==  2 && return 3.6953789770852417e-02 * 0.5
    i ==  3 && return 5.7134425426857240e-02 * 0.5
    i ==  4 && return 7.6100113628379332e-02 * 0.5
    i ==  5 && return 9.3444423456033876e-02 * 0.5
    i ==  6 && return 1.0879729916714846e-01 * 0.5
    i ==  7 && return 1.2183141605372846e-01 * 0.5
    i ==  8 && return 1.3226893863333750e-01 * 0.5
    i ==  9 && return 1.3988739479107321e-01 * 0.5
    i == 10 && return 1.4452440398997013e-01 * 0.5
    i == 11 && return 1.4608113364969041e-01 * 0.5
    i == 12 && return 1.4452440398997013e-01 * 0.5
    i == 13 && return 1.3988739479107321e-01 * 0.5
    i == 14 && return 1.3226893863333750e-01 * 0.5
    i == 15 && return 1.2183141605372846e-01 * 0.5
    i == 16 && return 1.0879729916714846e-01 * 0.5
    i == 17 && return 9.3444423456033876e-02 * 0.5
    i == 18 && return 7.6100113628379332e-02 * 0.5
    i == 19 && return 5.7134425426857240e-02 * 0.5
    i == 20 && return 3.6953789770852417e-02 * 0.5
    i == 21 && return 1.6017228257773849e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{21} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{21}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9375217062038945e-01,) .* 0.5
    i ==  2 && return (-9.6722683856630631e-01,) .* 0.5
    i ==  3 && return (-9.2009933415040079e-01,) .* 0.5
    i ==  4 && return (-8.5336336458331730e-01,) .* 0.5
    i ==  5 && return (-7.6843996347567789e-01,) .* 0.5
    i ==  6 && return (-6.6713880419741234e-01,) .* 0.5
    i ==  7 && return (-5.5161883588721983e-01,) .* 0.5
    i ==  8 && return (-4.2434212020743878e-01,) .* 0.5
    i ==  9 && return (-2.8802131680240112e-01,) .* 0.5
    i == 10 && return (-1.4556185416089509e-01,) .* 0.5
    i == 11 && return ( 0.0000000000000000e+00,) .* 0.5
    i == 12 && return ( 1.4556185416089509e-01,) .* 0.5
    i == 13 && return ( 2.8802131680240112e-01,) .* 0.5
    i == 14 && return ( 4.2434212020743878e-01,) .* 0.5
    i == 15 && return ( 5.5161883588721983e-01,) .* 0.5
    i == 16 && return ( 6.6713880419741234e-01,) .* 0.5
    i == 17 && return ( 7.6843996347567789e-01,) .* 0.5
    i == 18 && return ( 8.5336336458331730e-01,) .* 0.5
    i == 19 && return ( 9.2009933415040079e-01,) .* 0.5
    i == 20 && return ( 9.6722683856630631e-01,) .* 0.5
    i == 21 && return ( 9.9375217062038945e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{21} in AbstractEdge"))
end

# --------------------------- Strength φ = 22 ---------------------------
# Number of points: 22

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{22}) = 22

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{22}, i::Integer)::Float64
    i ==  1 && return 1.4627995298272575e-02 * 0.5
    i ==  2 && return 3.3774901584814200e-02 * 0.5
    i ==  3 && return 5.2293335152683272e-02 * 0.5
    i ==  4 && return 6.9796468424520447e-02 * 0.5
    i ==  5 && return 8.5941606217067840e-02 * 0.5
    i ==  6 && return 1.0041414444288101e-01 * 0.5
    i ==  7 && return 1.1293229608053923e-01 * 0.5
    i ==  8 && return 1.2325237681051242e-01 * 0.5
    i ==  9 && return 1.3117350478706241e-01 * 0.5
    i == 10 && return 1.3654149834601514e-01 * 0.5
    i == 11 && return 1.3925187285563204e-01 * 0.5
    i == 12 && return 1.3925187285563204e-01 * 0.5
    i == 13 && return 1.3654149834601514e-01 * 0.5
    i == 14 && return 1.3117350478706241e-01 * 0.5
    i == 15 && return 1.2325237681051242e-01 * 0.5
    i == 16 && return 1.1293229608053923e-01 * 0.5
    i == 17 && return 1.0041414444288101e-01 * 0.5
    i == 18 && return 8.5941606217067840e-02 * 0.5
    i == 19 && return 6.9796468424520447e-02 * 0.5
    i == 20 && return 5.2293335152683272e-02 * 0.5
    i == 21 && return 3.3774901584814200e-02 * 0.5
    i == 22 && return 1.4627995298272575e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{22} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{22}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9429458548239935e-01,) .* 0.5
    i ==  2 && return (-9.7006049783542869e-01,) .* 0.5
    i ==  3 && return (-9.2695677218717398e-01,) .* 0.5
    i ==  4 && return (-8.6581257772030018e-01,) .* 0.5
    i ==  5 && return (-7.8781680597920811e-01,) .* 0.5
    i ==  6 && return (-6.9448726318668275e-01,) .* 0.5
    i ==  7 && return (-5.8764040350691160e-01,) .* 0.5
    i ==  8 && return (-4.6935583798675701e-01,) .* 0.5
    i ==  9 && return (-3.4193582089208424e-01,) .* 0.5
    i == 10 && return (-2.0786042668822127e-01,) .* 0.5
    i == 11 && return (-6.9739273319722225e-02,) .* 0.5
    i == 12 && return ( 6.9739273319722225e-02,) .* 0.5
    i == 13 && return ( 2.0786042668822127e-01,) .* 0.5
    i == 14 && return ( 3.4193582089208424e-01,) .* 0.5
    i == 15 && return ( 4.6935583798675701e-01,) .* 0.5
    i == 16 && return ( 5.8764040350691160e-01,) .* 0.5
    i == 17 && return ( 6.9448726318668275e-01,) .* 0.5
    i == 18 && return ( 7.8781680597920811e-01,) .* 0.5
    i == 19 && return ( 8.6581257772030018e-01,) .* 0.5
    i == 20 && return ( 9.2695677218717398e-01,) .* 0.5
    i == 21 && return ( 9.7006049783542869e-01,) .* 0.5
    i == 22 && return ( 9.9429458548239935e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{22} in AbstractEdge"))
end

# --------------------------- Strength φ = 23 ---------------------------
# Number of points: 23

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{23}) = 23

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{23}, i::Integer)::Float64
    i ==  1 && return 1.3411859487141724e-02 * 0.5
    i ==  2 && return 3.0988005856979448e-02 * 0.5
    i ==  3 && return 4.8037671731084627e-02 * 0.5
    i ==  4 && return 6.4232421408525822e-02 * 0.5
    i ==  5 && return 7.9281411776718935e-02 * 0.5
    i ==  6 && return 9.2915766060035182e-02 * 0.5
    i ==  7 && return 1.0489209146454143e-01 * 0.5
    i ==  8 && return 1.1499664022241136e-01 * 0.5
    i ==  9 && return 1.2304908430672962e-01 * 0.5
    i == 10 && return 1.2890572218808219e-01 * 0.5
    i == 11 && return 1.3246203940469664e-01 * 0.5
    i == 12 && return 1.3365457218610619e-01 * 0.5
    i == 13 && return 1.3246203940469664e-01 * 0.5
    i == 14 && return 1.2890572218808219e-01 * 0.5
    i == 15 && return 1.2304908430672962e-01 * 0.5
    i == 16 && return 1.1499664022241136e-01 * 0.5
    i == 17 && return 1.0489209146454143e-01 * 0.5
    i == 18 && return 9.2915766060035182e-02 * 0.5
    i == 19 && return 7.9281411776718935e-02 * 0.5
    i == 20 && return 6.4232421408525822e-02 * 0.5
    i == 21 && return 4.8037671731084627e-02 * 0.5
    i == 22 && return 3.0988005856979448e-02 * 0.5
    i == 23 && return 1.3411859487141724e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{23} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{23}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9476933499755216e-01,) .* 0.5
    i ==  2 && return (-9.7254247121811521e-01,) .* 0.5
    i ==  3 && return (-9.3297108682601615e-01,) .* 0.5
    i ==  4 && return (-8.7675235827044162e-01,) .* 0.5
    i ==  5 && return (-8.0488840161883990e-01,) .* 0.5
    i ==  6 && return (-7.1866136313195017e-01,) .* 0.5
    i ==  7 && return (-6.1960987576364612e-01,) .* 0.5
    i ==  8 && return (-5.0950147784600752e-01,) .* 0.5
    i ==  9 && return (-3.9030103803029081e-01,) .* 0.5
    i == 10 && return (-2.6413568097034495e-01,) .* 0.5
    i == 11 && return (-1.3325682429846611e-01,) .* 0.5
    i == 12 && return ( 0.0000000000000000e+00,) .* 0.5
    i == 13 && return ( 1.3325682429846611e-01,) .* 0.5
    i == 14 && return ( 2.6413568097034495e-01,) .* 0.5
    i == 15 && return ( 3.9030103803029081e-01,) .* 0.5
    i == 16 && return ( 5.0950147784600752e-01,) .* 0.5
    i == 17 && return ( 6.1960987576364612e-01,) .* 0.5
    i == 18 && return ( 7.1866136313195017e-01,) .* 0.5
    i == 19 && return ( 8.0488840161883990e-01,) .* 0.5
    i == 20 && return ( 8.7675235827044162e-01,) .* 0.5
    i == 21 && return ( 9.3297108682601615e-01,) .* 0.5
    i == 22 && return ( 9.7254247121811521e-01,) .* 0.5
    i == 23 && return ( 9.9476933499755216e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{23} in AbstractEdge"))
end

# --------------------------- Strength φ = 24 ---------------------------
# Number of points: 24

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{24}) = 24

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{24}, i::Integer)::Float64
    i ==  1 && return 1.2341229799987641e-02 * 0.5
    i ==  2 && return 2.8531388628933584e-02 * 0.5
    i ==  3 && return 4.4277438817419801e-02 * 0.5
    i ==  4 && return 5.9298584915436693e-02 * 0.5
    i ==  5 && return 7.3346481411080272e-02 * 0.5
    i ==  6 && return 8.6190161531953191e-02 * 0.5
    i ==  7 && return 9.7618652104113898e-02 * 0.5
    i ==  8 && return 1.0744427011596554e-01 * 0.5
    i ==  9 && return 1.1550566805372557e-01 * 0.5
    i == 10 && return 1.2167047292780338e-01 * 0.5
    i == 11 && return 1.2583745634682833e-01 * 0.5
    i == 12 && return 1.2793819534675219e-01 * 0.5
    i == 13 && return 1.2793819534675219e-01 * 0.5
    i == 14 && return 1.2583745634682833e-01 * 0.5
    i == 15 && return 1.2167047292780338e-01 * 0.5
    i == 16 && return 1.1550566805372557e-01 * 0.5
    i == 17 && return 1.0744427011596554e-01 * 0.5
    i == 18 && return 9.7618652104113898e-02 * 0.5
    i == 19 && return 8.6190161531953191e-02 * 0.5
    i == 20 && return 7.3346481411080272e-02 * 0.5
    i == 21 && return 5.9298584915436693e-02 * 0.5
    i == 22 && return 4.4277438817419801e-02 * 0.5
    i == 23 && return 2.8531388628933584e-02 * 0.5
    i == 24 && return 1.2341229799987641e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{24} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{24}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9518721999702142e-01,) .* 0.5
    i ==  2 && return (-9.7472855597130947e-01,) .* 0.5
    i ==  3 && return (-9.3827455200273280e-01,) .* 0.5
    i ==  4 && return (-8.8641552700440107e-01,) .* 0.5
    i ==  5 && return (-8.2000198597390295e-01,) .* 0.5
    i ==  6 && return (-7.4012419157855436e-01,) .* 0.5
    i ==  7 && return (-6.4809365193697555e-01,) .* 0.5
    i ==  8 && return (-5.4542147138883956e-01,) .* 0.5
    i ==  9 && return (-4.3379350762604513e-01,) .* 0.5
    i == 10 && return (-3.1504267969616334e-01,) .* 0.5
    i == 11 && return (-1.9111886747361631e-01,) .* 0.5
    i == 12 && return (-6.4056892862605630e-02,) .* 0.5
    i == 13 && return ( 6.4056892862605630e-02,) .* 0.5
    i == 14 && return ( 1.9111886747361631e-01,) .* 0.5
    i == 15 && return ( 3.1504267969616334e-01,) .* 0.5
    i == 16 && return ( 4.3379350762604513e-01,) .* 0.5
    i == 17 && return ( 5.4542147138883956e-01,) .* 0.5
    i == 18 && return ( 6.4809365193697555e-01,) .* 0.5
    i == 19 && return ( 7.4012419157855436e-01,) .* 0.5
    i == 20 && return ( 8.2000198597390295e-01,) .* 0.5
    i == 21 && return ( 8.8641552700440107e-01,) .* 0.5
    i == 22 && return ( 9.3827455200273280e-01,) .* 0.5
    i == 23 && return ( 9.7472855597130947e-01,) .* 0.5
    i == 24 && return ( 9.9518721999702142e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{24} in AbstractEdge"))
end

# --------------------------- Strength φ = 25 ---------------------------
# Number of points: 25

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{25}) = 25

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{25}, i::Integer)::Float64
    i ==  1 && return 1.1393798501026250e-02 * 0.5
    i ==  2 && return 2.6354986615031935e-02 * 0.5
    i ==  3 && return 4.0939156701306247e-02 * 0.5
    i ==  4 && return 5.4904695975835305e-02 * 0.5
    i ==  5 && return 6.8038333812356952e-02 * 0.5
    i ==  6 && return 8.0140700335001133e-02 * 0.5
    i ==  7 && return 9.1028261982963626e-02 * 0.5
    i ==  8 && return 1.0053594906705073e-01 * 0.5
    i ==  9 && return 1.0851962447426364e-01 * 0.5
    i == 10 && return 1.1485825914571161e-01 * 0.5
    i == 11 && return 1.1945576353578476e-01 * 0.5
    i == 12 && return 1.2224244299031012e-01 * 0.5
    i == 13 && return 1.2317605372671545e-01 * 0.5
    i == 14 && return 1.2224244299031012e-01 * 0.5
    i == 15 && return 1.1945576353578476e-01 * 0.5
    i == 16 && return 1.1485825914571161e-01 * 0.5
    i == 17 && return 1.0851962447426364e-01 * 0.5
    i == 18 && return 1.0053594906705073e-01 * 0.5
    i == 19 && return 9.1028261982963626e-02 * 0.5
    i == 20 && return 8.0140700335001133e-02 * 0.5
    i == 21 && return 6.8038333812356952e-02 * 0.5
    i == 22 && return 5.4904695975835305e-02 * 0.5
    i == 23 && return 4.0939156701306247e-02 * 0.5
    i == 24 && return 2.6354986615031935e-02 * 0.5
    i == 25 && return 1.1393798501026250e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{25} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{25}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9555696979049813e-01,) .* 0.5
    i ==  2 && return (-9.7666392145951753e-01,) .* 0.5
    i ==  3 && return (-9.4297457122897432e-01,) .* 0.5
    i ==  4 && return (-8.9499199787827532e-01,) .* 0.5
    i ==  5 && return (-8.3344262876083397e-01,) .* 0.5
    i ==  6 && return (-7.5925926303735769e-01,) .* 0.5
    i ==  7 && return (-6.7356636847346840e-01,) .* 0.5
    i ==  8 && return (-5.7766293024122295e-01,) .* 0.5
    i ==  9 && return (-4.7300273144571497e-01,) .* 0.5
    i == 10 && return (-3.6117230580938786e-01,) .* 0.5
    i == 11 && return (-2.4386688372098844e-01,) .* 0.5
    i == 12 && return (-1.2286469261071040e-01,) .* 0.5
    i == 13 && return ( 0.0000000000000000e+00,) .* 0.5
    i == 14 && return ( 1.2286469261071040e-01,) .* 0.5
    i == 15 && return ( 2.4386688372098844e-01,) .* 0.5
    i == 16 && return ( 3.6117230580938786e-01,) .* 0.5
    i == 17 && return ( 4.7300273144571497e-01,) .* 0.5
    i == 18 && return ( 5.7766293024122295e-01,) .* 0.5
    i == 19 && return ( 6.7356636847346840e-01,) .* 0.5
    i == 20 && return ( 7.5925926303735769e-01,) .* 0.5
    i == 21 && return ( 8.3344262876083397e-01,) .* 0.5
    i == 22 && return ( 8.9499199787827532e-01,) .* 0.5
    i == 23 && return ( 9.4297457122897432e-01,) .* 0.5
    i == 24 && return ( 9.7666392145951753e-01,) .* 0.5
    i == 25 && return ( 9.9555696979049813e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{25} in AbstractEdge"))
end

# --------------------------- Strength φ = 26 ---------------------------
# Number of points: 26

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{26}) = 26

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{26}, i::Integer)::Float64
    i ==  1 && return 1.0551372617343093e-02 * 0.5
    i ==  2 && return 2.4417851092631858e-02 * 0.5
    i ==  3 && return 3.7962383294362800e-02 * 0.5
    i ==  4 && return 5.0975825297147767e-02 * 0.5
    i ==  5 && return 6.3274046329574812e-02 * 0.5
    i ==  6 && return 7.4684149765659749e-02 * 0.5
    i ==  7 && return 8.5045894313485151e-02 * 0.5
    i ==  8 && return 9.4213800355914118e-02 * 0.5
    i ==  9 && return 1.0205916109442549e-01 * 0.5
    i == 10 && return 1.0847184052857659e-01 * 0.5
    i == 11 && return 1.1336181654631960e-01 * 0.5
    i == 12 && return 1.1666044348529653e-01 * 0.5
    i == 13 && return 1.1832141527926224e-01 * 0.5
    i == 14 && return 1.1832141527926224e-01 * 0.5
    i == 15 && return 1.1666044348529653e-01 * 0.5
    i == 16 && return 1.1336181654631960e-01 * 0.5
    i == 17 && return 1.0847184052857659e-01 * 0.5
    i == 18 && return 1.0205916109442549e-01 * 0.5
    i == 19 && return 9.4213800355914118e-02 * 0.5
    i == 20 && return 8.5045894313485151e-02 * 0.5
    i == 21 && return 7.4684149765659749e-02 * 0.5
    i == 22 && return 6.3274046329574812e-02 * 0.5
    i == 23 && return 5.0975825297147767e-02 * 0.5
    i == 24 && return 3.7962383294362800e-02 * 0.5
    i == 25 && return 2.4417851092631858e-02 * 0.5
    i == 26 && return 1.0551372617343093e-02 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{26} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{26}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9588570114561692e-01,) .* 0.5
    i ==  2 && return (-9.7838544595647103e-01,) .* 0.5
    i ==  3 && return (-9.4715906666171423e-01,) .* 0.5
    i ==  4 && return (-9.0263786198430707e-01,) .* 0.5
    i ==  5 && return (-8.4544594278849805e-01,) .* 0.5
    i ==  6 && return (-7.7638594882067891e-01,) .* 0.5
    i ==  7 && return (-6.9642726041995728e-01,) .* 0.5
    i ==  8 && return (-6.0669229301761807e-01,) .* 0.5
    i ==  9 && return (-5.0844071482450570e-01,) .* 0.5
    i == 10 && return (-4.0305175512348629e-01,) .* 0.5
    i == 11 && return (-2.9200483948595690e-01,) .* 0.5
    i == 12 && return (-1.7685882035689018e-01,) .* 0.5
    i == 13 && return (-5.9230093429313208e-02,) .* 0.5
    i == 14 && return ( 5.9230093429313208e-02,) .* 0.5
    i == 15 && return ( 1.7685882035689018e-01,) .* 0.5
    i == 16 && return ( 2.9200483948595690e-01,) .* 0.5
    i == 17 && return ( 4.0305175512348629e-01,) .* 0.5
    i == 18 && return ( 5.0844071482450570e-01,) .* 0.5
    i == 19 && return ( 6.0669229301761807e-01,) .* 0.5
    i == 20 && return ( 6.9642726041995728e-01,) .* 0.5
    i == 21 && return ( 7.7638594882067891e-01,) .* 0.5
    i == 22 && return ( 8.4544594278849805e-01,) .* 0.5
    i == 23 && return ( 9.0263786198430707e-01,) .* 0.5
    i == 24 && return ( 9.4715906666171423e-01,) .* 0.5
    i == 25 && return ( 9.7838544595647103e-01,) .* 0.5
    i == 26 && return ( 9.9588570114561692e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{26} in AbstractEdge"))
end

# --------------------------- Strength φ = 27 ---------------------------
# Number of points: 27

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{27}) = 27

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{27}, i::Integer)::Float64
    i ==  1 && return 9.7989960512948962e-03 * 0.5
    i ==  2 && return 2.2686231596180564e-02 * 0.5
    i ==  3 && return 3.5297053757419726e-02 * 0.5
    i ==  4 && return 4.7449412520615145e-02 * 0.5
    i ==  5 && return 5.8983536859833562e-02 * 0.5
    i ==  6 && return 6.9748823766245527e-02 * 0.5
    i ==  7 && return 7.9604867773057766e-02 * 0.5
    i ==  8 && return 8.8423158543756986e-02 * 0.5
    i ==  9 && return 9.6088727370028590e-02 * 0.5
    i == 10 && return 1.0250163781774571e-01 * 0.5
    i == 11 && return 1.0757828578853322e-01 * 0.5
    i == 12 && return 1.1125248835684519e-01 * 0.5
    i == 13 && return 1.1347634610896513e-01 * 0.5
    i == 14 && return 1.1422086737895699e-01 * 0.5
    i == 15 && return 1.1347634610896513e-01 * 0.5
    i == 16 && return 1.1125248835684519e-01 * 0.5
    i == 17 && return 1.0757828578853322e-01 * 0.5
    i == 18 && return 1.0250163781774571e-01 * 0.5
    i == 19 && return 9.6088727370028590e-02 * 0.5
    i == 20 && return 8.8423158543756986e-02 * 0.5
    i == 21 && return 7.9604867773057766e-02 * 0.5
    i == 22 && return 6.9748823766245527e-02 * 0.5
    i == 23 && return 5.8983536859833562e-02 * 0.5
    i == 24 && return 4.7449412520615145e-02 * 0.5
    i == 25 && return 3.5297053757419726e-02 * 0.5
    i == 26 && return 2.2686231596180564e-02 * 0.5
    i == 27 && return 9.7989960512948962e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{27} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{27}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9617926288898861e-01,) .* 0.5
    i ==  2 && return (-9.7992347596150120e-01,) .* 0.5
    i ==  3 && return (-9.5090055781470495e-01,) .* 0.5
    i ==  4 && return (-9.0948232067749113e-01,) .* 0.5
    i ==  5 && return (-8.5620790801829449e-01,) .* 0.5
    i ==  6 && return (-7.9177163907050818e-01,) .* 0.5
    i ==  7 && return (-7.1701347373942370e-01,) .* 0.5
    i ==  8 && return (-6.3290797194649517e-01,) .* 0.5
    i ==  9 && return (-5.4055156457945686e-01,) .* 0.5
    i == 10 && return (-4.4114825175002687e-01,) .* 0.5
    i == 11 && return (-3.3599390363850890e-01,) .* 0.5
    i == 12 && return (-2.2645936543953685e-01,) .* 0.5
    i == 13 && return (-1.1397258560952997e-01,) .* 0.5
    i == 14 && return ( 0.0000000000000000e+00,) .* 0.5
    i == 15 && return ( 1.1397258560952997e-01,) .* 0.5
    i == 16 && return ( 2.2645936543953685e-01,) .* 0.5
    i == 17 && return ( 3.3599390363850890e-01,) .* 0.5
    i == 18 && return ( 4.4114825175002687e-01,) .* 0.5
    i == 19 && return ( 5.4055156457945686e-01,) .* 0.5
    i == 20 && return ( 6.3290797194649517e-01,) .* 0.5
    i == 21 && return ( 7.1701347373942370e-01,) .* 0.5
    i == 22 && return ( 7.9177163907050818e-01,) .* 0.5
    i == 23 && return ( 8.5620790801829449e-01,) .* 0.5
    i == 24 && return ( 9.0948232067749113e-01,) .* 0.5
    i == 25 && return ( 9.5090055781470495e-01,) .* 0.5
    i == 26 && return ( 9.7992347596150120e-01,) .* 0.5
    i == 27 && return ( 9.9617926288898861e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{27} in AbstractEdge"))
end

# --------------------------- Strength φ = 28 ---------------------------
# Number of points: 28

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{28}) = 28

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{28}, i::Integer)::Float64
    i ==  1 && return 9.1242825930945397e-03 * 0.5
    i ==  2 && return 2.1132112592771153e-02 * 0.5
    i ==  3 && return 3.2901427782304482e-02 * 0.5
    i ==  4 && return 4.4272934759004137e-02 * 0.5
    i ==  5 && return 5.5107345675716658e-02 * 0.5
    i ==  6 && return 6.5272923966999602e-02 * 0.5
    i ==  7 && return 7.4646214234568756e-02 * 0.5
    i ==  8 && return 8.3113417228901212e-02 * 0.5
    i ==  9 && return 9.0571744393032838e-02 * 0.5
    i == 10 && return 9.6930657997929950e-02 * 0.5
    i == 11 && return 1.0211296757806081e-01 * 0.5
    i == 12 && return 1.0605576592284642e-01 * 0.5
    i == 13 && return 1.0871119225829408e-01 * 0.5
    i == 14 && return 1.1004701301647522e-01 * 0.5
    i == 15 && return 1.1004701301647522e-01 * 0.5
    i == 16 && return 1.0871119225829408e-01 * 0.5
    i == 17 && return 1.0605576592284642e-01 * 0.5
    i == 18 && return 1.0211296757806081e-01 * 0.5
    i == 19 && return 9.6930657997929950e-02 * 0.5
    i == 20 && return 9.0571744393032838e-02 * 0.5
    i == 21 && return 8.3113417228901212e-02 * 0.5
    i == 22 && return 7.4646214234568756e-02 * 0.5
    i == 23 && return 6.5272923966999602e-02 * 0.5
    i == 24 && return 5.5107345675716658e-02 * 0.5
    i == 25 && return 4.4272934759004137e-02 * 0.5
    i == 26 && return 3.2901427782304482e-02 * 0.5
    i == 27 && return 2.1132112592771153e-02 * 0.5
    i == 28 && return 9.1242825930945397e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{28} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{28}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9644249757395442e-01,) .* 0.5
    i ==  2 && return (-9.8130316537087270e-01,) .* 0.5
    i ==  3 && return (-9.5425928062893817e-01,) .* 0.5
    i ==  4 && return (-9.1563302639213207e-01,) .* 0.5
    i ==  5 && return (-8.6589252257439508e-01,) .* 0.5
    i ==  6 && return (-8.0564137091717913e-01,) .* 0.5
    i ==  7 && return (-7.3561087801363179e-01,) .* 0.5
    i ==  8 && return (-6.5665109403886490e-01,) .* 0.5
    i ==  9 && return (-5.6972047181140173e-01,) .* 0.5
    i == 10 && return (-4.7587422495511827e-01,) .* 0.5
    i == 11 && return (-3.7625151608907870e-01,) .* 0.5
    i == 12 && return (-2.7206162763517810e-01,) .* 0.5
    i == 13 && return (-1.6456928213338076e-01,) .* 0.5
    i == 14 && return (-5.5079289884034259e-02,) .* 0.5
    i == 15 && return ( 5.5079289884034259e-02,) .* 0.5
    i == 16 && return ( 1.6456928213338076e-01,) .* 0.5
    i == 17 && return ( 2.7206162763517810e-01,) .* 0.5
    i == 18 && return ( 3.7625151608907870e-01,) .* 0.5
    i == 19 && return ( 4.7587422495511827e-01,) .* 0.5
    i == 20 && return ( 5.6972047181140173e-01,) .* 0.5
    i == 21 && return ( 6.5665109403886490e-01,) .* 0.5
    i == 22 && return ( 7.3561087801363179e-01,) .* 0.5
    i == 23 && return ( 8.0564137091717913e-01,) .* 0.5
    i == 24 && return ( 8.6589252257439508e-01,) .* 0.5
    i == 25 && return ( 9.1563302639213207e-01,) .* 0.5
    i == 26 && return ( 9.5425928062893817e-01,) .* 0.5
    i == 27 && return ( 9.8130316537087270e-01,) .* 0.5
    i == 28 && return ( 9.9644249757395442e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{28} in AbstractEdge"))
end

# --------------------------- Strength φ = 29 ---------------------------
# Number of points: 29

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{29}) = 29

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{29}, i::Integer)::Float64
    i ==  1 && return 8.5169038787464933e-03 * 0.5
    i ==  2 && return 1.9732085056122738e-02 * 0.5
    i ==  3 && return 3.0740492202093652e-02 * 0.5
    i ==  4 && return 4.1402062518682906e-02 * 0.5
    i ==  5 && return 5.1594826902497885e-02 * 0.5
    i ==  6 && return 6.1203090657079094e-02 * 0.5
    i ==  7 && return 7.0117933255051210e-02 * 0.5
    i ==  8 && return 7.8238327135763744e-02 * 0.5
    i ==  9 && return 8.5472257366172560e-02 * 0.5
    i == 10 && return 9.1737757139258733e-02 * 0.5
    i == 11 && return 9.6963834094408674e-02 * 0.5
    i == 12 && return 1.0109127375991502e-01 * 0.5
    i == 13 && return 1.0407331007772935e-01 * 0.5
    i == 14 && return 1.0587615509732100e-01 * 0.5
    i == 15 && return 1.0647938171831425e-01 * 0.5
    i == 16 && return 1.0587615509732100e-01 * 0.5
    i == 17 && return 1.0407331007772935e-01 * 0.5
    i == 18 && return 1.0109127375991502e-01 * 0.5
    i == 19 && return 9.6963834094408674e-02 * 0.5
    i == 20 && return 9.1737757139258733e-02 * 0.5
    i == 21 && return 8.5472257366172560e-02 * 0.5
    i == 22 && return 7.8238327135763744e-02 * 0.5
    i == 23 && return 7.0117933255051210e-02 * 0.5
    i == 24 && return 6.1203090657079094e-02 * 0.5
    i == 25 && return 5.1594826902497885e-02 * 0.5
    i == 26 && return 4.1402062518682906e-02 * 0.5
    i == 27 && return 3.0740492202093652e-02 * 0.5
    i == 28 && return 1.9732085056122738e-02 * 0.5
    i == 29 && return 8.5169038787464933e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{29} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{29}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9667944226059657e-01,) .* 0.5
    i ==  2 && return (-9.8254550526141315e-01,) .* 0.5
    i ==  3 && return (-9.5728559577808769e-01,) .* 0.5
    i ==  4 && return (-9.2118023295305873e-01,) .* 0.5
    i ==  5 && return (-8.7463780492010279e-01,) .* 0.5
    i ==  6 && return (-8.1818548761525245e-01,) .* 0.5
    i ==  7 && return (-7.5246285173447713e-01,) .* 0.5
    i ==  8 && return (-6.7821453760268646e-01,) .* 0.5
    i ==  9 && return (-5.9628179713822782e-01,) .* 0.5
    i == 10 && return (-5.0759295512422764e-01,) .* 0.5
    i == 11 && return (-4.1315288817400869e-01,) .* 0.5
    i == 12 && return (-3.1403163786763993e-01,) .* 0.5
    i == 13 && return (-2.1135228616600107e-01,) .* 0.5
    i == 14 && return (-1.0627823013267922e-01,) .* 0.5
    i == 15 && return ( 0.0000000000000000e+00,) .* 0.5
    i == 16 && return ( 1.0627823013267922e-01,) .* 0.5
    i == 17 && return ( 2.1135228616600107e-01,) .* 0.5
    i == 18 && return ( 3.1403163786763993e-01,) .* 0.5
    i == 19 && return ( 4.1315288817400869e-01,) .* 0.5
    i == 20 && return ( 5.0759295512422764e-01,) .* 0.5
    i == 21 && return ( 5.9628179713822782e-01,) .* 0.5
    i == 22 && return ( 6.7821453760268646e-01,) .* 0.5
    i == 23 && return ( 7.5246285173447713e-01,) .* 0.5
    i == 24 && return ( 8.1818548761525245e-01,) .* 0.5
    i == 25 && return ( 8.7463780492010279e-01,) .* 0.5
    i == 26 && return ( 9.2118023295305873e-01,) .* 0.5
    i == 27 && return ( 9.5728559577808769e-01,) .* 0.5
    i == 28 && return ( 9.8254550526141315e-01,) .* 0.5
    i == 29 && return ( 9.9667944226059657e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{29} in AbstractEdge"))
end

# --------------------------- Strength φ = 30 ---------------------------
# Number of points: 30

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{30}) = 30

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{30}, i::Integer)::Float64
    i ==  1 && return 7.9681924961665183e-03 * 0.5
    i ==  2 && return 1.8466468311090854e-02 * 0.5
    i ==  3 && return 2.8784707883323272e-02 * 0.5
    i ==  4 && return 3.8799192569627022e-02 * 0.5
    i ==  5 && return 4.8402672830594073e-02 * 0.5
    i ==  6 && return 5.7493156217619086e-02 * 0.5
    i ==  7 && return 6.5974229882180518e-02 * 0.5
    i ==  8 && return 7.3755974737705149e-02 * 0.5
    i ==  9 && return 8.0755895229420227e-02 * 0.5
    i == 10 && return 8.6899787201082962e-02 * 0.5
    i == 11 && return 9.2122522237786136e-02 * 0.5
    i == 12 && return 9.6368737174644198e-02 * 0.5
    i == 13 && return 9.9593420586795281e-02 * 0.5
    i == 14 && return 1.0176238974840554e-01 * 0.5
    i == 15 && return 1.0285265289355890e-01 * 0.5
    i == 16 && return 1.0285265289355890e-01 * 0.5
    i == 17 && return 1.0176238974840554e-01 * 0.5
    i == 18 && return 9.9593420586795281e-02 * 0.5
    i == 19 && return 9.6368737174644198e-02 * 0.5
    i == 20 && return 9.2122522237786136e-02 * 0.5
    i == 21 && return 8.6899787201082962e-02 * 0.5
    i == 22 && return 8.0755895229420227e-02 * 0.5
    i == 23 && return 7.3755974737705149e-02 * 0.5
    i == 24 && return 6.5974229882180518e-02 * 0.5
    i == 25 && return 5.7493156217619086e-02 * 0.5
    i == 26 && return 4.8402672830594073e-02 * 0.5
    i == 27 && return 3.8799192569627022e-02 * 0.5
    i == 28 && return 2.8784707883323272e-02 * 0.5
    i == 29 && return 1.8466468311090854e-02 * 0.5
    i == 30 && return 7.9681924961665183e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{30} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{30}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9689348407464962e-01,) .* 0.5
    i ==  2 && return (-9.8366812327974718e-01,) .* 0.5
    i ==  3 && return (-9.6002186496830755e-01,) .* 0.5
    i ==  4 && return (-9.2620004742927431e-01,) .* 0.5
    i ==  5 && return (-8.8256053579205263e-01,) .* 0.5
    i ==  6 && return (-8.2956576238276836e-01,) .* 0.5
    i ==  7 && return (-7.6777743210482619e-01,) .* 0.5
    i ==  8 && return (-6.9785049479331585e-01,) .* 0.5
    i ==  9 && return (-6.2052618298924289e-01,) .* 0.5
    i == 10 && return (-5.3662414814201986e-01,) .* 0.5
    i == 11 && return (-4.4703376953808915e-01,) .* 0.5
    i == 12 && return (-3.5270472553087812e-01,) .* 0.5
    i == 13 && return (-2.5463692616788985e-01,) .* 0.5
    i == 14 && return (-1.5386991360858354e-01,) .* 0.5
    i == 15 && return (-5.1471842555317691e-02,) .* 0.5
    i == 16 && return ( 5.1471842555317691e-02,) .* 0.5
    i == 17 && return ( 1.5386991360858354e-01,) .* 0.5
    i == 18 && return ( 2.5463692616788985e-01,) .* 0.5
    i == 19 && return ( 3.5270472553087812e-01,) .* 0.5
    i == 20 && return ( 4.4703376953808915e-01,) .* 0.5
    i == 21 && return ( 5.3662414814201986e-01,) .* 0.5
    i == 22 && return ( 6.2052618298924289e-01,) .* 0.5
    i == 23 && return ( 6.9785049479331585e-01,) .* 0.5
    i == 24 && return ( 7.6777743210482619e-01,) .* 0.5
    i == 25 && return ( 8.2956576238276836e-01,) .* 0.5
    i == 26 && return ( 8.8256053579205263e-01,) .* 0.5
    i == 27 && return ( 9.2620004742927431e-01,) .* 0.5
    i == 28 && return ( 9.6002186496830755e-01,) .* 0.5
    i == 29 && return ( 9.8366812327974718e-01,) .* 0.5
    i == 30 && return ( 9.9689348407464962e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{30} in AbstractEdge"))
end

# --------------------------- Strength φ = 31 ---------------------------
# Number of points: 31

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{31}) = 31

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{31}, i::Integer)::Float64
    i ==  1 && return 7.4708315792487460e-03 * 0.5
    i ==  2 && return 1.7318620790310619e-02 * 0.5
    i ==  3 && return 2.7009019184979371e-02 * 0.5
    i ==  4 && return 3.6432273912385509e-02 * 0.5
    i ==  5 && return 4.5493707527201006e-02 * 0.5
    i ==  6 && return 5.4103082424916842e-02 * 0.5
    i ==  7 && return 6.2174786561028456e-02 * 0.5
    i ==  8 && return 6.9628583235410380e-02 * 0.5
    i ==  9 && return 7.6390386598776602e-02 * 0.5
    i == 10 && return 8.2392991761589235e-02 * 0.5
    i == 11 && return 8.7576740608477782e-02 * 0.5
    i == 12 && return 9.1890113893641504e-02 * 0.5
    i == 13 && return 9.5290242912319523e-02 * 0.5
    i == 14 && return 9.7743335386328775e-02 * 0.5
    i == 15 && return 9.9225011226672266e-02 * 0.5
    i == 16 && return 9.9720544793426444e-02 * 0.5
    i == 17 && return 9.9225011226672266e-02 * 0.5
    i == 18 && return 9.7743335386328775e-02 * 0.5
    i == 19 && return 9.5290242912319523e-02 * 0.5
    i == 20 && return 9.1890113893641504e-02 * 0.5
    i == 21 && return 8.7576740608477782e-02 * 0.5
    i == 22 && return 8.2392991761589235e-02 * 0.5
    i == 23 && return 7.6390386598776602e-02 * 0.5
    i == 24 && return 6.9628583235410380e-02 * 0.5
    i == 25 && return 6.2174786561028456e-02 * 0.5
    i == 26 && return 5.4103082424916842e-02 * 0.5
    i == 27 && return 4.5493707527201006e-02 * 0.5
    i == 28 && return 3.6432273912385509e-02 * 0.5
    i == 29 && return 2.7009019184979371e-02 * 0.5
    i == 30 && return 1.7318620790310619e-02 * 0.5
    i == 31 && return 7.4708315792487460e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{31} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{31}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9708748181947715e-01,) .* 0.5
    i ==  2 && return (-9.8468590966515246e-01,) .* 0.5
    i ==  3 && return (-9.6250392509294969e-01,) .* 0.5
    i ==  4 && return (-9.3075699789664812e-01,) .* 0.5
    i ==  5 && return (-8.8976002994827108e-01,) .* 0.5
    i ==  6 && return (-8.3992032014626739e-01,) .* 0.5
    i ==  7 && return (-7.8173314841662500e-01,) .* 0.5
    i ==  8 && return (-7.1577678458685323e-01,) .* 0.5
    i ==  9 && return (-6.4270672292426034e-01,) .* 0.5
    i == 10 && return (-5.6324916140714931e-01,) .* 0.5
    i == 11 && return (-4.7819378204490248e-01,) .* 0.5
    i == 12 && return (-3.8838590160823294e-01,) .* 0.5
    i == 13 && return (-2.9471806998170164e-01,) .* 0.5
    i == 14 && return (-1.9812119933557062e-01,) .* 0.5
    i == 15 && return (-9.9555312152341521e-02,) .* 0.5
    i == 16 && return ( 0.0000000000000000e+00,) .* 0.5
    i == 17 && return ( 9.9555312152341521e-02,) .* 0.5
    i == 18 && return ( 1.9812119933557062e-01,) .* 0.5
    i == 19 && return ( 2.9471806998170164e-01,) .* 0.5
    i == 20 && return ( 3.8838590160823294e-01,) .* 0.5
    i == 21 && return ( 4.7819378204490248e-01,) .* 0.5
    i == 22 && return ( 5.6324916140714931e-01,) .* 0.5
    i == 23 && return ( 6.4270672292426034e-01,) .* 0.5
    i == 24 && return ( 7.1577678458685323e-01,) .* 0.5
    i == 25 && return ( 7.8173314841662500e-01,) .* 0.5
    i == 26 && return ( 8.3992032014626739e-01,) .* 0.5
    i == 27 && return ( 8.8976002994827108e-01,) .* 0.5
    i == 28 && return ( 9.3075699789664812e-01,) .* 0.5
    i == 29 && return ( 9.6250392509294969e-01,) .* 0.5
    i == 30 && return ( 9.8468590966515246e-01,) .* 0.5
    i == 31 && return ( 9.9708748181947715e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{31} in AbstractEdge"))
end

# --------------------------- Strength φ = 32 ---------------------------
# Number of points: 32

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{32}) = 32

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{32}, i::Integer)::Float64
    i ==  1 && return 7.0186100094700920e-03 * 0.5
    i ==  2 && return 1.6274394730905643e-02 * 0.5
    i ==  3 && return 2.5392065309262066e-02 * 0.5
    i ==  4 && return 3.4273862913021411e-02 * 0.5
    i ==  5 && return 4.2835898022226718e-02 * 0.5
    i ==  6 && return 5.0998059262376147e-02 * 0.5
    i ==  7 && return 5.8684093478535482e-02 * 0.5
    i ==  8 && return 6.5822222776361794e-02 * 0.5
    i ==  9 && return 7.2345794108848560e-02 * 0.5
    i == 10 && return 7.8193895787070325e-02 * 0.5
    i == 11 && return 8.3311924226946776e-02 * 0.5
    i == 12 && return 8.7652093004403769e-02 * 0.5
    i == 13 && return 9.1173878695763835e-02 * 0.5
    i == 14 && return 9.3844399080804650e-02 * 0.5
    i == 15 && return 9.5638720079274889e-02 * 0.5
    i == 16 && return 9.6540088514727784e-02 * 0.5
    i == 17 && return 9.6540088514727784e-02 * 0.5
    i == 18 && return 9.5638720079274889e-02 * 0.5
    i == 19 && return 9.3844399080804650e-02 * 0.5
    i == 20 && return 9.1173878695763835e-02 * 0.5
    i == 21 && return 8.7652093004403769e-02 * 0.5
    i == 22 && return 8.3311924226946776e-02 * 0.5
    i == 23 && return 7.8193895787070325e-02 * 0.5
    i == 24 && return 7.2345794108848560e-02 * 0.5
    i == 25 && return 6.5822222776361794e-02 * 0.5
    i == 26 && return 5.8684093478535482e-02 * 0.5
    i == 27 && return 5.0998059262376147e-02 * 0.5
    i == 28 && return 4.2835898022226718e-02 * 0.5
    i == 29 && return 3.4273862913021411e-02 * 0.5
    i == 30 && return 2.5392065309262066e-02 * 0.5
    i == 31 && return 1.6274394730905643e-02 * 0.5
    i == 32 && return 7.0186100094700920e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{32} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{32}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9726386184948157e-01,) .* 0.5
    i ==  2 && return (-9.8561151154526838e-01,) .* 0.5
    i ==  3 && return (-9.6476225558750639e-01,) .* 0.5
    i ==  4 && return (-9.3490607593773967e-01,) .* 0.5
    i ==  5 && return (-8.9632115576605209e-01,) .* 0.5
    i ==  6 && return (-8.4936761373256997e-01,) .* 0.5
    i ==  7 && return (-7.9448379596794239e-01,) .* 0.5
    i ==  8 && return (-7.3218211874028971e-01,) .* 0.5
    i ==  9 && return (-6.6304426693021523e-01,) .* 0.5
    i == 10 && return (-5.8771575724076230e-01,) .* 0.5
    i == 11 && return (-5.0689990893222936e-01,) .* 0.5
    i == 12 && return (-4.2135127613063533e-01,) .* 0.5
    i == 13 && return (-3.3186860228212767e-01,) .* 0.5
    i == 14 && return (-2.3928736225213706e-01,) .* 0.5
    i == 15 && return (-1.4447196158279649e-01,) .* 0.5
    i == 16 && return (-4.8307665687738324e-02,) .* 0.5
    i == 17 && return ( 4.8307665687738324e-02,) .* 0.5
    i == 18 && return ( 1.4447196158279649e-01,) .* 0.5
    i == 19 && return ( 2.3928736225213706e-01,) .* 0.5
    i == 20 && return ( 3.3186860228212767e-01,) .* 0.5
    i == 21 && return ( 4.2135127613063533e-01,) .* 0.5
    i == 22 && return ( 5.0689990893222936e-01,) .* 0.5
    i == 23 && return ( 5.8771575724076230e-01,) .* 0.5
    i == 24 && return ( 6.6304426693021523e-01,) .* 0.5
    i == 25 && return ( 7.3218211874028971e-01,) .* 0.5
    i == 26 && return ( 7.9448379596794239e-01,) .* 0.5
    i == 27 && return ( 8.4936761373256997e-01,) .* 0.5
    i == 28 && return ( 8.9632115576605209e-01,) .* 0.5
    i == 29 && return ( 9.3490607593773967e-01,) .* 0.5
    i == 30 && return ( 9.6476225558750639e-01,) .* 0.5
    i == 31 && return ( 9.8561151154526838e-01,) .* 0.5
    i == 32 && return ( 9.9726386184948157e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{32} in AbstractEdge"))
end

# --------------------------- Strength φ = 33 ---------------------------
# Number of points: 33

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{33}) = 33

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{33}, i::Integer)::Float64
    i ==  1 && return 6.6062278475873554e-03 * 0.5
    i ==  2 && return 1.5321701512934717e-02 * 0.5
    i ==  3 && return 2.3915548101749520e-02 * 0.5
    i ==  4 && return 3.2300358632328940e-02 * 0.5
    i ==  5 && return 4.0401541331669628e-02 * 0.5
    i ==  6 && return 4.8147742818711793e-02 * 0.5
    i ==  7 && return 5.5470846631663566e-02 * 0.5
    i ==  8 && return 6.2306482530317502e-02 * 0.5
    i ==  9 && return 6.8594572818656704e-02 * 0.5
    i == 10 && return 7.4279854843954107e-02 * 0.5
    i == 11 && return 7.9312364794886750e-02 * 0.5
    i == 12 && return 8.3647876067038676e-02 * 0.5
    i == 13 && return 8.7248287618844358e-02 * 0.5
    i == 14 && return 9.0081958660638534e-02 * 0.5
    i == 15 && return 9.2123986643316863e-02 * 0.5
    i == 16 && return 9.3356426065596049e-02 * 0.5
    i == 17 && return 9.3768446160209989e-02 * 0.5
    i == 18 && return 9.3356426065596049e-02 * 0.5
    i == 19 && return 9.2123986643316863e-02 * 0.5
    i == 20 && return 9.0081958660638534e-02 * 0.5
    i == 21 && return 8.7248287618844358e-02 * 0.5
    i == 22 && return 8.3647876067038676e-02 * 0.5
    i == 23 && return 7.9312364794886750e-02 * 0.5
    i == 24 && return 7.4279854843954107e-02 * 0.5
    i == 25 && return 6.8594572818656704e-02 * 0.5
    i == 26 && return 6.2306482530317502e-02 * 0.5
    i == 27 && return 5.5470846631663566e-02 * 0.5
    i == 28 && return 4.8147742818711793e-02 * 0.5
    i == 29 && return 4.0401541331669628e-02 * 0.5
    i == 30 && return 3.2300358632328940e-02 * 0.5
    i == 31 && return 2.3915548101749520e-02 * 0.5
    i == 32 && return 1.5321701512934717e-02 * 0.5
    i == 33 && return 6.6062278475873554e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{33} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{33}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9742469424645519e-01,) .* 0.5
    i ==  2 && return (-9.8645572623064248e-01,) .* 0.5
    i ==  3 && return (-9.6682290968999274e-01,) .* 0.5
    i ==  4 && return (-9.3869437261116839e-01,) .* 0.5
    i ==  5 && return (-9.0231676774343361e-01,) .* 0.5
    i ==  6 && return (-8.5800965267650409e-01,) .* 0.5
    i ==  7 && return (-8.0616235627416655e-01,) .* 0.5
    i ==  8 && return (-7.4723049644956219e-01,) .* 0.5
    i ==  9 && return (-6.8173195996974278e-01,) .* 0.5
    i == 10 && return (-6.1024234583637904e-01,) .* 0.5
    i == 11 && return (-5.3338990478634762e-01,) .* 0.5
    i == 12 && return (-4.5185001727245072e-01,) .* 0.5
    i == 13 && return (-3.6633925774807335e-01,) .* 0.5
    i == 14 && return (-2.7760909715249704e-01,) .* 0.5
    i == 15 && return (-1.8643929882799157e-01,) .* 0.5
    i == 16 && return (-9.3631065854733381e-02,) .* 0.5
    i == 17 && return ( 0.0000000000000000e+00,) .* 0.5
    i == 18 && return ( 9.3631065854733381e-02,) .* 0.5
    i == 19 && return ( 1.8643929882799157e-01,) .* 0.5
    i == 20 && return ( 2.7760909715249704e-01,) .* 0.5
    i == 21 && return ( 3.6633925774807335e-01,) .* 0.5
    i == 22 && return ( 4.5185001727245072e-01,) .* 0.5
    i == 23 && return ( 5.3338990478634762e-01,) .* 0.5
    i == 24 && return ( 6.1024234583637904e-01,) .* 0.5
    i == 25 && return ( 6.8173195996974278e-01,) .* 0.5
    i == 26 && return ( 7.4723049644956219e-01,) .* 0.5
    i == 27 && return ( 8.0616235627416655e-01,) .* 0.5
    i == 28 && return ( 8.5800965267650409e-01,) .* 0.5
    i == 29 && return ( 9.0231676774343361e-01,) .* 0.5
    i == 30 && return ( 9.3869437261116839e-01,) .* 0.5
    i == 31 && return ( 9.6682290968999274e-01,) .* 0.5
    i == 32 && return ( 9.8645572623064248e-01,) .* 0.5
    i == 33 && return ( 9.9742469424645519e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{33} in AbstractEdge"))
end

# --------------------------- Strength φ = 34 ---------------------------
# Number of points: 34

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{34}) = 34

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{34}, i::Integer)::Float64
    i ==  1 && return 6.2291405559085937e-03 * 0.5
    i ==  2 && return 1.4450162748595181e-02 * 0.5
    i ==  3 && return 2.2563721985495035e-02 * 0.5
    i ==  4 && return 3.0491380638446110e-02 * 0.5
    i ==  5 && return 3.8166593796387559e-02 * 0.5
    i ==  6 && return 4.5525611523353250e-02 * 0.5
    i ==  7 && return 5.2507414572678122e-02 * 0.5
    i ==  8 && return 5.9054135827524522e-02 * 0.5
    i ==  9 && return 6.5111521554076457e-02 * 0.5
    i == 10 && return 7.0629375814255782e-02 * 0.5
    i == 11 && return 7.5561974660031908e-02 * 0.5
    i == 12 && return 7.9868444339771832e-02 * 0.5
    i == 13 && return 8.3513099699845661e-02 * 0.5
    i == 14 && return 8.6465739747035794e-02 * 0.5
    i == 15 && return 8.8701897835693780e-02 * 0.5
    i == 16 && return 9.0203044370640736e-02 * 0.5
    i == 17 && return 9.0956740330259842e-02 * 0.5
    i == 18 && return 9.0956740330259842e-02 * 0.5
    i == 19 && return 9.0203044370640736e-02 * 0.5
    i == 20 && return 8.8701897835693780e-02 * 0.5
    i == 21 && return 8.6465739747035794e-02 * 0.5
    i == 22 && return 8.3513099699845661e-02 * 0.5
    i == 23 && return 7.9868444339771832e-02 * 0.5
    i == 24 && return 7.5561974660031908e-02 * 0.5
    i == 25 && return 7.0629375814255782e-02 * 0.5
    i == 26 && return 6.5111521554076457e-02 * 0.5
    i == 27 && return 5.9054135827524522e-02 * 0.5
    i == 28 && return 5.2507414572678122e-02 * 0.5
    i == 29 && return 4.5525611523353250e-02 * 0.5
    i == 30 && return 3.8166593796387559e-02 * 0.5
    i == 31 && return 3.0491380638446110e-02 * 0.5
    i == 32 && return 2.2563721985495035e-02 * 0.5
    i == 33 && return 1.4450162748595181e-02 * 0.5
    i == 34 && return 6.2291405559085937e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{34} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{34}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9757175379084195e-01,) .* 0.5
    i ==  2 && return (-9.8722781640630952e-01,) .* 0.5
    i ==  3 && return (-9.6870826253334430e-01,) .* 0.5
    i ==  4 && return (-9.4216239740510710e-01,) .* 0.5
    i ==  5 && return (-9.0780967771832444e-01,) .* 0.5
    i ==  6 && return (-8.6593463833456452e-01,) .* 0.5
    i ==  7 && return (-8.1688422790093362e-01,) .* 0.5
    i ==  8 && return (-7.6106487662987299e-01,) .* 0.5
    i ==  9 && return (-6.9893911321626290e-01,) .* 0.5
    i == 10 && return (-6.3102172708052851e-01,) .* 0.5
    i == 11 && return (-5.5787550066974667e-01,) .* 0.5
    i == 12 && return (-4.8010654519032703e-01,) .* 0.5
    i == 13 && return (-3.9835927775864594e-01,) .* 0.5
    i == 14 && return (-3.1331108133946323e-01,) .* 0.5
    i == 15 && return (-2.2566669161644948e-01,) .* 0.5
    i == 16 && return (-1.3615235725918298e-01,) .* 0.5
    i == 17 && return (-4.5509821953102547e-02,) .* 0.5
    i == 18 && return ( 4.5509821953102547e-02,) .* 0.5
    i == 19 && return ( 1.3615235725918298e-01,) .* 0.5
    i == 20 && return ( 2.2566669161644948e-01,) .* 0.5
    i == 21 && return ( 3.1331108133946323e-01,) .* 0.5
    i == 22 && return ( 3.9835927775864594e-01,) .* 0.5
    i == 23 && return ( 4.8010654519032703e-01,) .* 0.5
    i == 24 && return ( 5.5787550066974667e-01,) .* 0.5
    i == 25 && return ( 6.3102172708052851e-01,) .* 0.5
    i == 26 && return ( 6.9893911321626290e-01,) .* 0.5
    i == 27 && return ( 7.6106487662987299e-01,) .* 0.5
    i == 28 && return ( 8.1688422790093362e-01,) .* 0.5
    i == 29 && return ( 8.6593463833456452e-01,) .* 0.5
    i == 30 && return ( 9.0780967771832444e-01,) .* 0.5
    i == 31 && return ( 9.4216239740510710e-01,) .* 0.5
    i == 32 && return ( 9.6870826253334430e-01,) .* 0.5
    i == 33 && return ( 9.8722781640630952e-01,) .* 0.5
    i == 34 && return ( 9.9757175379084195e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{34} in AbstractEdge"))
end

# --------------------------- Strength φ = 35 ---------------------------
# Number of points: 35

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{35}) = 35

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{35}, i::Integer)::Float64
    i ==  1 && return 5.8834334204429838e-03 * 0.5
    i ==  2 && return 1.3650828348361692e-02 * 0.5
    i ==  3 && return 2.1322979911483606e-02 * 0.5
    i ==  4 && return 2.8829260108894146e-02 * 0.5
    i ==  5 && return 3.6110115863463348e-02 * 0.5
    i ==  6 && return 4.3108422326170188e-02 * 0.5
    i ==  7 && return 4.9769370401353596e-02 * 0.5
    i ==  8 && return 5.6040816212370108e-02 * 0.5
    i ==  9 && return 6.1873671966080221e-02 * 0.5
    i == 10 && return 6.7222285269086926e-02 * 0.5
    i == 11 && return 7.2044794772560122e-02 * 0.5
    i == 12 && return 7.6303457155441998e-02 * 0.5
    i == 13 && return 7.9964942242324269e-02 * 0.5
    i == 14 && return 8.3000593728856611e-02 * 0.5
    i == 15 && return 8.5386653392099152e-02 * 0.5
    i == 16 && return 8.7104446997183491e-02 * 0.5
    i == 17 && return 8.8140530430275504e-02 * 0.5
    i == 18 && return 8.8486794907104288e-02 * 0.5
    i == 19 && return 8.8140530430275504e-02 * 0.5
    i == 20 && return 8.7104446997183491e-02 * 0.5
    i == 21 && return 8.5386653392099152e-02 * 0.5
    i == 22 && return 8.3000593728856611e-02 * 0.5
    i == 23 && return 7.9964942242324269e-02 * 0.5
    i == 24 && return 7.6303457155441998e-02 * 0.5
    i == 25 && return 7.2044794772560122e-02 * 0.5
    i == 26 && return 6.7222285269086926e-02 * 0.5
    i == 27 && return 6.1873671966080221e-02 * 0.5
    i == 28 && return 5.6040816212370108e-02 * 0.5
    i == 29 && return 4.9769370401353596e-02 * 0.5
    i == 30 && return 4.3108422326170188e-02 * 0.5
    i == 31 && return 3.6110115863463348e-02 * 0.5
    i == 32 && return 2.8829260108894146e-02 * 0.5
    i == 33 && return 2.1322979911483606e-02 * 0.5
    i == 34 && return 1.3650828348361692e-02 * 0.5
    i == 35 && return 5.8834334204429838e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{35} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{35}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9770656909960032e-01,) .* 0.5
    i ==  2 && return (-9.8793576444385156e-01,) .* 0.5
    i ==  3 && return (-9.7043761603922984e-01,) .* 0.5
    i ==  4 && return (-9.4534514820782733e-01,) .* 0.5
    i ==  5 && return (-9.1285426135931758e-01,) .* 0.5
    i ==  6 && return (-8.7321912502522236e-01,) .* 0.5
    i ==  7 && return (-8.2674989909222540e-01,) .* 0.5
    i ==  8 && return (-7.7381025228691258e-01,) .* 0.5
    i ==  9 && return (-7.1481450155662885e-01,) .* 0.5
    i == 10 && return (-6.5022436466589040e-01,) .* 0.5
    i == 11 && return (-5.8054534474976449e-01,) .* 0.5
    i == 12 && return (-5.0632277324148867e-01,) .* 0.5
    i == 13 && return (-4.2813754151781425e-01,) .* 0.5
    i == 14 && return (-3.4660155443081392e-01,) .* 0.5
    i == 15 && return (-2.6235294120929603e-01,) .* 0.5
    i == 16 && return (-1.7605106116598956e-01,) .* 0.5
    i == 17 && return (-8.8371343275659264e-02,) .* 0.5
    i == 18 && return ( 0.0000000000000000e+00,) .* 0.5
    i == 19 && return ( 8.8371343275659264e-02,) .* 0.5
    i == 20 && return ( 1.7605106116598956e-01,) .* 0.5
    i == 21 && return ( 2.6235294120929603e-01,) .* 0.5
    i == 22 && return ( 3.4660155443081392e-01,) .* 0.5
    i == 23 && return ( 4.2813754151781425e-01,) .* 0.5
    i == 24 && return ( 5.0632277324148867e-01,) .* 0.5
    i == 25 && return ( 5.8054534474976449e-01,) .* 0.5
    i == 26 && return ( 6.5022436466589040e-01,) .* 0.5
    i == 27 && return ( 7.1481450155662885e-01,) .* 0.5
    i == 28 && return ( 7.7381025228691258e-01,) .* 0.5
    i == 29 && return ( 8.2674989909222540e-01,) .* 0.5
    i == 30 && return ( 8.7321912502522236e-01,) .* 0.5
    i == 31 && return ( 9.1285426135931758e-01,) .* 0.5
    i == 32 && return ( 9.4534514820782733e-01,) .* 0.5
    i == 33 && return ( 9.7043761603922984e-01,) .* 0.5
    i == 34 && return ( 9.8793576444385156e-01,) .* 0.5
    i == 35 && return ( 9.9770656909960032e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{35} in AbstractEdge"))
end

# --------------------------- Strength φ = 36 ---------------------------
# Number of points: 36

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{36}) = 36

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{36}, i::Integer)::Float64
    i ==  1 && return 5.5657196642451617e-03 * 0.5
    i ==  2 && return 1.2915947284065629e-02 * 0.5
    i ==  3 && return 2.0181515297735594e-02 * 0.5
    i ==  4 && return 2.7298621498568747e-02 * 0.5
    i ==  5 && return 3.4213810770307281e-02 * 0.5
    i ==  6 && return 4.0875750923644975e-02 * 0.5
    i ==  7 && return 4.7235083490265950e-02 * 0.5
    i ==  8 && return 5.3244713977759872e-02 * 0.5
    i ==  9 && return 5.8860144245324937e-02 * 0.5
    i == 10 && return 6.4039797355015471e-02 * 0.5
    i == 11 && return 6.8745323835736463e-02 * 0.5
    i == 12 && return 7.2941885005653059e-02 * 0.5
    i == 13 && return 7.6598410645870640e-02 * 0.5
    i == 14 && return 7.9687828912071615e-02 * 0.5
    i == 15 && return 8.2187266704339651e-02 * 0.5
    i == 16 && return 8.4078218979661973e-02 * 0.5
    i == 17 && return 8.5346685739338665e-02 * 0.5
    i == 18 && return 8.5983275670394751e-02 * 0.5
    i == 19 && return 8.5983275670394751e-02 * 0.5
    i == 20 && return 8.5346685739338665e-02 * 0.5
    i == 21 && return 8.4078218979661973e-02 * 0.5
    i == 22 && return 8.2187266704339651e-02 * 0.5
    i == 23 && return 7.9687828912071615e-02 * 0.5
    i == 24 && return 7.6598410645870640e-02 * 0.5
    i == 25 && return 7.2941885005653059e-02 * 0.5
    i == 26 && return 6.8745323835736463e-02 * 0.5
    i == 27 && return 6.4039797355015471e-02 * 0.5
    i == 28 && return 5.8860144245324937e-02 * 0.5
    i == 29 && return 5.3244713977759872e-02 * 0.5
    i == 30 && return 4.7235083490265950e-02 * 0.5
    i == 31 && return 4.0875750923644975e-02 * 0.5
    i == 32 && return 3.4213810770307281e-02 * 0.5
    i == 33 && return 2.7298621498568747e-02 * 0.5
    i == 34 && return 2.0181515297735594e-02 * 0.5
    i == 35 && return 1.2915947284065629e-02 * 0.5
    i == 36 && return 5.5657196642451617e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{36} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{36}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9783046248408580e-01,) .* 0.5
    i ==  2 && return (-9.8858647890221218e-01,) .* 0.5
    i ==  3 && return (-9.7202769104969799e-01,) .* 0.5
    i ==  4 && return (-9.4827298439950758e-01,) .* 0.5
    i ==  5 && return (-9.1749777451565906e-01,) .* 0.5
    i ==  6 && return (-8.7992980089039707e-01,) .* 0.5
    i ==  7 && return (-8.3584716699247530e-01,) .* 0.5
    i ==  8 && return (-7.8557623013220657e-01,) .* 0.5
    i ==  9 && return (-7.2948917159355653e-01,) .* 0.5
    i == 10 && return (-6.6800123658552102e-01,) .* 0.5
    i == 11 && return (-6.0156765813598057e-01,) .* 0.5
    i == 12 && return (-5.3068028592624517e-01,) .* 0.5
    i == 13 && return (-4.5586394443342027e-01,) .* 0.5
    i == 14 && return (-3.7767254711968923e-01,) .* 0.5
    i == 15 && return (-2.9668499534402826e-01,) .* 0.5
    i == 16 && return (-2.1350089231686559e-01,) .* 0.5
    i == 17 && return (-1.2873610380938480e-01,) .* 0.5
    i == 18 && return (-4.3018198473708601e-02,) .* 0.5
    i == 19 && return ( 4.3018198473708601e-02,) .* 0.5
    i == 20 && return ( 1.2873610380938480e-01,) .* 0.5
    i == 21 && return ( 2.1350089231686559e-01,) .* 0.5
    i == 22 && return ( 2.9668499534402826e-01,) .* 0.5
    i == 23 && return ( 3.7767254711968923e-01,) .* 0.5
    i == 24 && return ( 4.5586394443342027e-01,) .* 0.5
    i == 25 && return ( 5.3068028592624517e-01,) .* 0.5
    i == 26 && return ( 6.0156765813598057e-01,) .* 0.5
    i == 27 && return ( 6.6800123658552102e-01,) .* 0.5
    i == 28 && return ( 7.2948917159355653e-01,) .* 0.5
    i == 29 && return ( 7.8557623013220657e-01,) .* 0.5
    i == 30 && return ( 8.3584716699247530e-01,) .* 0.5
    i == 31 && return ( 8.7992980089039707e-01,) .* 0.5
    i == 32 && return ( 9.1749777451565906e-01,) .* 0.5
    i == 33 && return ( 9.4827298439950758e-01,) .* 0.5
    i == 34 && return ( 9.7202769104969799e-01,) .* 0.5
    i == 35 && return ( 9.8858647890221218e-01,) .* 0.5
    i == 36 && return ( 9.9783046248408580e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{36} in AbstractEdge"))
end

# --------------------------- Strength φ = 37 ---------------------------
# Number of points: 37

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{37}) = 37

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{37}, i::Integer)::Float64
    i ==  1 && return 5.2730572794980656e-03 * 0.5
    i ==  2 && return 1.2238780100307479e-02 * 0.5
    i ==  3 && return 1.9129044489084007e-02 * 0.5
    i ==  4 && return 2.5886036990558990e-02 * 0.5
    i ==  5 && return 3.2461639847521484e-02 * 0.5
    i ==  6 && return 3.8809602501934423e-02 * 0.5
    i ==  7 && return 4.4885364662437192e-02 * 0.5
    i ==  8 && return 5.0646297654824556e-02 * 0.5
    i ==  9 && return 5.6051987998274877e-02 * 0.5
    i == 10 && return 6.1064516523226010e-02 * 0.5
    i == 11 && return 6.5648722872751170e-02 * 0.5
    i == 12 && return 6.9772451555700318e-02 * 0.5
    i == 13 && return 7.3406777248488181e-02 * 0.5
    i == 14 && return 7.6526207570529164e-02 * 0.5
    i == 15 && return 7.9108861837529409e-02 * 0.5
    i == 16 && return 8.1136624508465052e-02 * 0.5
    i == 17 && return 8.2595272236437256e-02 * 0.5
    i == 18 && return 8.3474573625862775e-02 * 0.5
    i == 19 && return 8.3768360993138904e-02 * 0.5
    i == 20 && return 8.3474573625862775e-02 * 0.5
    i == 21 && return 8.2595272236437256e-02 * 0.5
    i == 22 && return 8.1136624508465052e-02 * 0.5
    i == 23 && return 7.9108861837529409e-02 * 0.5
    i == 24 && return 7.6526207570529164e-02 * 0.5
    i == 25 && return 7.3406777248488181e-02 * 0.5
    i == 26 && return 6.9772451555700318e-02 * 0.5
    i == 27 && return 6.5648722872751170e-02 * 0.5
    i == 28 && return 6.1064516523226010e-02 * 0.5
    i == 29 && return 5.6051987998274877e-02 * 0.5
    i == 30 && return 5.0646297654824556e-02 * 0.5
    i == 31 && return 4.4885364662437192e-02 * 0.5
    i == 32 && return 3.8809602501934423e-02 * 0.5
    i == 33 && return 3.2461639847521484e-02 * 0.5
    i == 34 && return 2.5886036990558990e-02 * 0.5
    i == 35 && return 1.9129044489084007e-02 * 0.5
    i == 36 && return 1.2238780100307479e-02 * 0.5
    i == 37 && return 5.2730572794980656e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{37} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{37}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9794458247791362e-01,) .* 0.5
    i ==  2 && return (-9.8918596321431917e-01,) .* 0.5
    i ==  3 && return (-9.7349303005648569e-01,) .* 0.5
    i ==  4 && return (-9.5097234326209479e-01,) .* 0.5
    i ==  5 && return (-9.2178143741246377e-01,) .* 0.5
    i ==  6 && return (-8.8612496215548608e-01,) .* 0.5
    i ==  7 && return (-8.4425298734055598e-01,) .* 0.5
    i ==  8 && return (-7.9645920050990227e-01,) .* 0.5
    i ==  9 && return (-7.4307883398196528e-01,) .* 0.5
    i == 10 && return (-6.8448630913095931e-01,) .* 0.5
    i == 11 && return (-6.2109260840892444e-01,) .* 0.5
    i == 12 && return (-5.5334239186158174e-01,) .* 0.5
    i == 13 && return (-4.8171087780320554e-01,) .* 0.5
    i == 14 && return (-4.0670050931832613e-01,) .* 0.5
    i == 15 && return (-3.2883742988370701e-01,) .* 0.5
    i == 16 && return (-2.4866779279136575e-01,) .* 0.5
    i == 17 && return (-1.6675393023985197e-01,) .* 0.5
    i == 18 && return (-8.3670408954769890e-02,) .* 0.5
    i == 19 && return ( 0.0000000000000000e+00,) .* 0.5
    i == 20 && return ( 8.3670408954769890e-02,) .* 0.5
    i == 21 && return ( 1.6675393023985197e-01,) .* 0.5
    i == 22 && return ( 2.4866779279136575e-01,) .* 0.5
    i == 23 && return ( 3.2883742988370701e-01,) .* 0.5
    i == 24 && return ( 4.0670050931832613e-01,) .* 0.5
    i == 25 && return ( 4.8171087780320554e-01,) .* 0.5
    i == 26 && return ( 5.5334239186158174e-01,) .* 0.5
    i == 27 && return ( 6.2109260840892444e-01,) .* 0.5
    i == 28 && return ( 6.8448630913095931e-01,) .* 0.5
    i == 29 && return ( 7.4307883398196528e-01,) .* 0.5
    i == 30 && return ( 7.9645920050990227e-01,) .* 0.5
    i == 31 && return ( 8.4425298734055598e-01,) .* 0.5
    i == 32 && return ( 8.8612496215548608e-01,) .* 0.5
    i == 33 && return ( 9.2178143741246377e-01,) .* 0.5
    i == 34 && return ( 9.5097234326209479e-01,) .* 0.5
    i == 35 && return ( 9.7349303005648569e-01,) .* 0.5
    i == 36 && return ( 9.8918596321431917e-01,) .* 0.5
    i == 37 && return ( 9.9794458247791362e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{37} in AbstractEdge"))
end

# --------------------------- Strength φ = 38 ---------------------------
# Number of points: 38

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{38}) = 38

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{38}, i::Integer)::Float64
    i ==  1 && return 5.0028807496393353e-03 * 0.5
    i ==  2 && return 1.1613444716468704e-02 * 0.5
    i ==  3 && return 1.8156577709613233e-02 * 0.5
    i ==  4 && return 2.4579739738232350e-02 * 0.5
    i ==  5 && return 3.0839500545175026e-02 * 0.5
    i ==  6 && return 3.6894081594024748e-02 * 0.5
    i ==  7 && return 4.2703158504674432e-02 * 0.5
    i ==  8 && return 4.8228061860758682e-02 * 0.5
    i ==  9 && return 5.3432019910332307e-02 * 0.5
    i == 10 && return 5.8280399146997175e-02 * 0.5
    i == 11 && return 6.2740933392133075e-02 * 0.5
    i == 12 && return 6.6783937979140368e-02 * 0.5
    i == 13 && return 7.0382507066898983e-02 * 0.5
    i == 14 && return 7.3512692584743397e-02 * 0.5
    i == 15 && return 7.6153663548446424e-02 * 0.5
    i == 16 && return 7.8287844658210898e-02 * 0.5
    i == 17 && return 7.9901033243527833e-02 * 0.5
    i == 18 && return 8.0982493770597061e-02 * 0.5
    i == 19 && return 8.1525029280385755e-02 * 0.5
    i == 20 && return 8.1525029280385755e-02 * 0.5
    i == 21 && return 8.0982493770597061e-02 * 0.5
    i == 22 && return 7.9901033243527833e-02 * 0.5
    i == 23 && return 7.8287844658210898e-02 * 0.5
    i == 24 && return 7.6153663548446424e-02 * 0.5
    i == 25 && return 7.3512692584743397e-02 * 0.5
    i == 26 && return 7.0382507066898983e-02 * 0.5
    i == 27 && return 6.6783937979140368e-02 * 0.5
    i == 28 && return 6.2740933392133075e-02 * 0.5
    i == 29 && return 5.8280399146997175e-02 * 0.5
    i == 30 && return 5.3432019910332307e-02 * 0.5
    i == 31 && return 4.8228061860758682e-02 * 0.5
    i == 32 && return 4.2703158504674432e-02 * 0.5
    i == 33 && return 3.6894081594024748e-02 * 0.5
    i == 34 && return 3.0839500545175026e-02 * 0.5
    i == 35 && return 2.4579739738232350e-02 * 0.5
    i == 36 && return 1.8156577709613233e-02 * 0.5
    i == 37 && return 1.1613444716468704e-02 * 0.5
    i == 38 && return 5.0028807496393353e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{38} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{38}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9804993053568758e-01,) .* 0.5
    i ==  2 && return (-9.8973945426638554e-01,) .* 0.5
    i ==  3 && return (-9.7484632859015352e-01,) .* 0.5
    i ==  4 && return (-9.5346633093352962e-01,) .* 0.5
    i ==  5 && return (-9.2574133204858444e-01,) .* 0.5
    i ==  6 && return (-8.9185573900463222e-01,) .* 0.5
    i ==  7 && return (-8.5203502193236214e-01,) .* 0.5
    i ==  8 && return (-8.0654416760531678e-01,) .* 0.5
    i ==  9 && return (-7.5568590375397071e-01,) .* 0.5
    i == 10 && return (-6.9979868037918436e-01,) .* 0.5
    i == 11 && return (-6.3925441582968168e-01,) .* 0.5
    i == 12 && return (-5.7445602104780713e-01,) .* 0.5
    i == 13 && return (-5.0583471792793111e-01,) .* 0.5
    i == 14 && return (-4.3384716943237650e-01,) .* 0.5
    i == 15 && return (-3.5897244047943500e-01,) .* 0.5
    i == 16 && return (-2.8170880979016527e-01,) .* 0.5
    i == 17 && return (-2.0257045389211673e-01,) .* 0.5
    i == 18 && return (-1.2208402533786741e-01,) .* 0.5
    i == 19 && return (-4.0785147904578253e-02,) .* 0.5
    i == 20 && return ( 4.0785147904578253e-02,) .* 0.5
    i == 21 && return ( 1.2208402533786741e-01,) .* 0.5
    i == 22 && return ( 2.0257045389211673e-01,) .* 0.5
    i == 23 && return ( 2.8170880979016527e-01,) .* 0.5
    i == 24 && return ( 3.5897244047943500e-01,) .* 0.5
    i == 25 && return ( 4.3384716943237650e-01,) .* 0.5
    i == 26 && return ( 5.0583471792793111e-01,) .* 0.5
    i == 27 && return ( 5.7445602104780713e-01,) .* 0.5
    i == 28 && return ( 6.3925441582968168e-01,) .* 0.5
    i == 29 && return ( 6.9979868037918436e-01,) .* 0.5
    i == 30 && return ( 7.5568590375397071e-01,) .* 0.5
    i == 31 && return ( 8.0654416760531678e-01,) .* 0.5
    i == 32 && return ( 8.5203502193236214e-01,) .* 0.5
    i == 33 && return ( 8.9185573900463222e-01,) .* 0.5
    i == 34 && return ( 9.2574133204858444e-01,) .* 0.5
    i == 35 && return ( 9.5346633093352962e-01,) .* 0.5
    i == 36 && return ( 9.7484632859015352e-01,) .* 0.5
    i == 37 && return ( 9.8973945426638554e-01,) .* 0.5
    i == 38 && return ( 9.9804993053568758e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{38} in AbstractEdge"))
end

# --------------------------- Strength φ = 39 ---------------------------
# Number of points: 39

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{39}) = 39

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{39}, i::Integer)::Float64
    i ==  1 && return 4.7529446916349700e-03 * 0.5
    i ==  2 && return 1.1034788939164517e-02 * 0.5
    i ==  3 && return 1.7256229093724935e-02 * 0.5
    i ==  4 && return 2.3369384832178208e-02 * 0.5
    i ==  5 && return 2.9334955983903403e-02 * 0.5
    i ==  6 && return 3.5115111498131360e-02 * 0.5
    i ==  7 && return 4.0673276847933787e-02 * 0.5
    i ==  8 && return 4.5974301108916656e-02 * 0.5
    i ==  9 && return 5.0984665292129465e-02 * 0.5
    i == 10 && return 5.5672690340916271e-02 * 0.5
    i == 11 && return 6.0008736088596117e-02 * 0.5
    i == 12 && return 6.3965388138682369e-02 * 0.5
    i == 13 && return 6.7517630966231298e-02 * 0.5
    i == 14 && return 7.0643005970608783e-02 * 0.5
    i == 15 && return 7.3321753414268623e-02 * 0.5
    i == 16 && return 7.5536937322835992e-02 * 0.5
    i == 17 && return 7.7274552544682018e-02 * 0.5
    i == 18 && return 7.8523613287371119e-02 * 0.5
    i == 19 && return 7.9276222568368443e-02 * 0.5
    i == 20 && return 7.9527622139442852e-02 * 0.5
    i == 21 && return 7.9276222568368443e-02 * 0.5
    i == 22 && return 7.8523613287371119e-02 * 0.5
    i == 23 && return 7.7274552544682018e-02 * 0.5
    i == 24 && return 7.5536937322835992e-02 * 0.5
    i == 25 && return 7.3321753414268623e-02 * 0.5
    i == 26 && return 7.0643005970608783e-02 * 0.5
    i == 27 && return 6.7517630966231298e-02 * 0.5
    i == 28 && return 6.3965388138682369e-02 * 0.5
    i == 29 && return 6.0008736088596117e-02 * 0.5
    i == 30 && return 5.5672690340916271e-02 * 0.5
    i == 31 && return 5.0984665292129465e-02 * 0.5
    i == 32 && return 4.5974301108916656e-02 * 0.5
    i == 33 && return 4.0673276847933787e-02 * 0.5
    i == 34 && return 3.5115111498131360e-02 * 0.5
    i == 35 && return 2.9334955983903403e-02 * 0.5
    i == 36 && return 2.3369384832178208e-02 * 0.5
    i == 37 && return 1.7256229093724935e-02 * 0.5
    i == 38 && return 1.1034788939164517e-02 * 0.5
    i == 39 && return 4.7529446916349700e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{39} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{39}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9814738306643291e-01,) .* 0.5
    i ==  2 && return (-9.9025153685468603e-01,) .* 0.5
    i ==  3 && return (-9.7609870933347109e-01,) .* 0.5
    i ==  4 && return (-9.5577521232465223e-01,) .* 0.5
    i ==  5 && return (-9.2940914848673828e-01,) .* 0.5
    i ==  6 && return (-8.9716711929299287e-01,) .* 0.5
    i ==  7 && return (-8.5925293799990610e-01,) .* 0.5
    i ==  8 && return (-8.1590629743014309e-01,) .* 0.5
    i ==  9 && return (-7.6740124293106349e-01,) .* 0.5
    i == 10 && return (-7.1404443589453470e-01,) .* 0.5
    i == 11 && return (-6.5617321343201096e-01,) .* 0.5
    i == 12 && return (-5.9415345495727800e-01,) .* 0.5
    i == 13 && return (-5.2837726866043744e-01,) .* 0.5
    i == 14 && return (-4.5926051230913606e-01,) .* 0.5
    i == 15 && return (-3.8724016397156147e-01,) .* 0.5
    i == 16 && return (-3.1277155924818589e-01,) .* 0.5
    i == 17 && return (-2.3632551246183575e-01,) .* 0.5
    i == 18 && return (-1.5838533999783780e-01,) .* 0.5
    i == 19 && return (-7.9443804608755469e-02,) .* 0.5
    i == 20 && return ( 0.0000000000000000e+00,) .* 0.5
    i == 21 && return ( 7.9443804608755469e-02,) .* 0.5
    i == 22 && return ( 1.5838533999783780e-01,) .* 0.5
    i == 23 && return ( 2.3632551246183575e-01,) .* 0.5
    i == 24 && return ( 3.1277155924818589e-01,) .* 0.5
    i == 25 && return ( 3.8724016397156147e-01,) .* 0.5
    i == 26 && return ( 4.5926051230913606e-01,) .* 0.5
    i == 27 && return ( 5.2837726866043744e-01,) .* 0.5
    i == 28 && return ( 5.9415345495727800e-01,) .* 0.5
    i == 29 && return ( 6.5617321343201096e-01,) .* 0.5
    i == 30 && return ( 7.1404443589453470e-01,) .* 0.5
    i == 31 && return ( 7.6740124293106349e-01,) .* 0.5
    i == 32 && return ( 8.1590629743014309e-01,) .* 0.5
    i == 33 && return ( 8.5925293799990610e-01,) .* 0.5
    i == 34 && return ( 8.9716711929299287e-01,) .* 0.5
    i == 35 && return ( 9.2940914848673828e-01,) .* 0.5
    i == 36 && return ( 9.5577521232465223e-01,) .* 0.5
    i == 37 && return ( 9.7609870933347109e-01,) .* 0.5
    i == 38 && return ( 9.9025153685468603e-01,) .* 0.5
    i == 39 && return ( 9.9814738306643291e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{39} in AbstractEdge"))
end

# --------------------------- Strength φ = 40 ---------------------------
# Number of points: 40

nquadpoints(::Type{<:AbstractEdge}, ::GaussQuadrature{40}) = 40

function quadweight(::Type{<:AbstractEdge}, ::GaussQuadrature{40}, i::Integer)::Float64
    i ==  1 && return 4.5212770985337521e-03 * 0.5
    i ==  2 && return 1.0498284531152837e-02 * 0.5
    i ==  3 && return 1.6421058381907903e-02 * 0.5
    i ==  4 && return 2.2245849194166965e-02 * 0.5
    i ==  5 && return 2.7937006980023403e-02 * 0.5
    i ==  6 && return 3.3460195282547851e-02 * 0.5
    i ==  7 && return 3.8782167974471933e-02 * 0.5
    i ==  8 && return 4.3870908185673262e-02 * 0.5
    i ==  9 && return 4.8695807635072211e-02 * 0.5
    i == 10 && return 5.3227846983936733e-02 * 0.5
    i == 11 && return 5.7439769099391517e-02 * 0.5
    i == 12 && return 6.1306242492928965e-02 * 0.5
    i == 13 && return 6.4804013456601098e-02 * 0.5
    i == 14 && return 6.7912045815233871e-02 * 0.5
    i == 15 && return 7.0611647391286794e-02 * 0.5
    i == 16 && return 7.2886582395804075e-02 * 0.5
    i == 17 && return 7.4723169057968261e-02 * 0.5
    i == 18 && return 7.6110361900626269e-02 * 0.5
    i == 19 && return 7.7039818164248028e-02 * 0.5
    i == 20 && return 7.7505947978424777e-02 * 0.5
    i == 21 && return 7.7505947978424777e-02 * 0.5
    i == 22 && return 7.7039818164248028e-02 * 0.5
    i == 23 && return 7.6110361900626269e-02 * 0.5
    i == 24 && return 7.4723169057968261e-02 * 0.5
    i == 25 && return 7.2886582395804075e-02 * 0.5
    i == 26 && return 7.0611647391286794e-02 * 0.5
    i == 27 && return 6.7912045815233871e-02 * 0.5
    i == 28 && return 6.4804013456601098e-02 * 0.5
    i == 29 && return 6.1306242492928965e-02 * 0.5
    i == 30 && return 5.7439769099391517e-02 * 0.5
    i == 31 && return 5.3227846983936733e-02 * 0.5
    i == 32 && return 4.8695807635072211e-02 * 0.5
    i == 33 && return 4.3870908185673262e-02 * 0.5
    i == 34 && return 3.8782167974471933e-02 * 0.5
    i == 35 && return 3.3460195282547851e-02 * 0.5
    i == 36 && return 2.7937006980023403e-02 * 0.5
    i == 37 && return 2.2245849194166965e-02 * 0.5
    i == 38 && return 1.6421058381907903e-02 * 0.5
    i == 39 && return 1.0498284531152837e-02 * 0.5
    i == 40 && return 4.5212770985337521e-03 * 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{40} in AbstractEdge"))
end

function quadpoint(::Type{<:AbstractEdge}, ::GaussQuadrature{40}, i::Integer)::SVector{1, Float64}
    i ==  1 && return (-9.9823770971055925e-01,) .* 0.5
    i ==  2 && return (-9.9072623869945697e-01,) .* 0.5
    i ==  3 && return (-9.7725994998377430e-01,) .* 0.5
    i ==  4 && return (-9.5791681921379168e-01,) .* 0.5
    i ==  5 && return (-9.3281280827867652e-01,) .* 0.5
    i ==  6 && return (-9.0209880696887434e-01,) .* 0.5
    i ==  7 && return (-8.6595950321225945e-01,) .* 0.5
    i ==  8 && return (-8.2461223083331170e-01,) .* 0.5
    i ==  9 && return (-7.7830565142651942e-01,) .* 0.5
    i == 10 && return (-7.2731825518992710e-01,) .* 0.5
    i == 11 && return (-6.7195668461417957e-01,) .* 0.5
    i == 12 && return (-6.1255388966798019e-01,) .* 0.5
    i == 13 && return (-5.4946712509512818e-01,) .* 0.5
    i == 14 && return (-4.8307580168617870e-01,) .* 0.5
    i == 15 && return (-4.1377920437160498e-01,) .* 0.5
    i == 16 && return (-3.4199409082575849e-01,) .* 0.5
    i == 17 && return (-2.6815218500725369e-01,) .* 0.5
    i == 18 && return (-1.9269758070137111e-01,) .* 0.5
    i == 19 && return (-1.1608407067525521e-01,) .* 0.5
    i == 20 && return (-3.8772417506050823e-02,) .* 0.5
    i == 21 && return ( 3.8772417506050823e-02,) .* 0.5
    i == 22 && return ( 1.1608407067525521e-01,) .* 0.5
    i == 23 && return ( 1.9269758070137111e-01,) .* 0.5
    i == 24 && return ( 2.6815218500725369e-01,) .* 0.5
    i == 25 && return ( 3.4199409082575849e-01,) .* 0.5
    i == 26 && return ( 4.1377920437160498e-01,) .* 0.5
    i == 27 && return ( 4.8307580168617870e-01,) .* 0.5
    i == 28 && return ( 5.4946712509512818e-01,) .* 0.5
    i == 29 && return ( 6.1255388966798019e-01,) .* 0.5
    i == 30 && return ( 6.7195668461417957e-01,) .* 0.5
    i == 31 && return ( 7.2731825518992710e-01,) .* 0.5
    i == 32 && return ( 7.7830565142651942e-01,) .* 0.5
    i == 33 && return ( 8.2461223083331170e-01,) .* 0.5
    i == 34 && return ( 8.6595950321225945e-01,) .* 0.5
    i == 35 && return ( 9.0209880696887434e-01,) .* 0.5
    i == 36 && return ( 9.3281280827867652e-01,) .* 0.5
    i == 37 && return ( 9.5791681921379168e-01,) .* 0.5
    i == 38 && return ( 9.7725994998377430e-01,) .* 0.5
    i == 39 && return ( 9.9072623869945697e-01,) .* 0.5
    i == 40 && return ( 9.9823770971055925e-01,) .* 0.5
    throw(ArgumentError("no quadrature point $i for GaussQuadrature{40} in AbstractEdge"))
end
