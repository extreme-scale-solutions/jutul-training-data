import Random

function VtoK_simple(v::T; α1=T(4), α2=T(-10)) where T
    exp(α1*v + α2)
end

function VtoK_stochastic(rng::Random.AbstractRNG, maxin::T, thisvalue::T; α1=T(4), α2=T(-9.8), α3=T(0.0531660692), maxoutinterval=Array{T}([1000, 10000])) where T
    randlow = -3*α3
    randhigh = 3*α3
    α1 = α1 + clamp(randn(rng, T) * α3, randlow, randhigh)
    α2 = α2 + clamp(randn(rng, T) * α3, randlow, randhigh)
    maxout = VtoK_simple(maxin; α1=α1, α2=α2)
    if maxout < maxoutinterval[1] || maxout > maxoutinterval[2]
        @warn("bad VtoK_stochastic parameters: input=$maxin output=$maxout α1=$α1 α2=$α2 α3=$α3")
    end
    VtoK_simple(thisvalue; α1=α1, α2=α2)
end

function VtoK_banded(v::T; α1=T(4), α2=T(4.25), α3=T(3000), α4=T(3.5), α5=T(0.01), α6=T(1.5)) where T
    if v >= α1
        K = exp(v-α2) * α3 / exp(α1-α2)
    elseif v >= α4
        K = α5 * exp(log(α3/α5)/(α1-α4)*(v-α4))
    else
        K = α5 * exp(v-α6)/exp(α4-α6)
    end
    return K
end

function downsample(v::Matrix{T}, factor::Int) where T
    return downsample(v, (factor, factor))
end

function downsample(v::Matrix{T}, factor::Tuple{Int, Int}) where T
    v_out_size = div.(size(v), factor)
    v_out = zeros(T, v_out_size)
    for i = 1:v_out_size[1]
        for j = 1:v_out_size[2]
            v_out[i,j] = mean(v[factor[1]*i-factor[1]+1:factor[1]*i, factor[2]*j-factor[2]+1:factor[2]*j])
        end
    end
    return v_out
end

function downsample(v::Array{T, 3}, factor::Tuple{Int, Int, Int}) where T
    v_out_size = div.(size(v), factor)
    v_out = zeros(T, v_out_size)
    for i = 1:v_out_size[1]
        for j = 1:v_out_size[2]
            for k = 1:v_out_size[3]
                v_out[i,j,k] = mean(v[factor[1]*i-factor[1]+1:factor[1]*i, factor[2]*j-factor[2]+1:factor[2]*j, factor[3]*k-factor[3]+1:factor[3]*k])
            end
        end
    end
    return v_out
end

### kozeny-carman, α from UK CCS appraisal project report

function ϕtoK(ϕ::T; α=T(48.63057324840764)) where T
    return ϕ^3 * (α / (T(1)-ϕ))^2
end

function Ktoϕ(K::T; α=T(48.63057324840764)) where T
    p = Polynomial([-K,2*K,-K, α^2])
    return minimum(real(roots(p)[findall(real(roots(p)).== roots(p))]))
end

### pad porosity to avoid overpressure
function padϕ(ϕ::Matrix{T}) where T
    return hcat(vcat(T(1e8)*ones(T, 1, size(ϕ,2)-1), ϕ[2:end-1,1:end-1], T(1e8)*ones(T, 1, size(ϕ,2)-1)),
    T(1e8)*ones(T, size(ϕ,1), 1))
end

function padϕ(ϕ::Array{T, 3}) where T
    rv = copy(ϕ)
    # the following comments assume x,y,z coordinates are longitude (west-east), latitude (north-south), depth.  but it's arbitrary.
    # west side
    rv[1,:,:]   .= 1e8
    # east side
    rv[end,:,:] .= 1e8
    # north side
    rv[:,1,:]   .= 1e8
    # south side
    rv[:,end,:] .= 1e8
    # under side
    rv[:,:,end] .= 1e8
    # leave the top side alone, it's supposed to be caprock
    rv
end
