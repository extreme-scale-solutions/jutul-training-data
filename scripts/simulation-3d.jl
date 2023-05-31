### Read the Compass slice from SLIM website and take out a 3D slice
## Please refer to the details in the MIT license of this repository and in the license of the Compass model
## author: Ziyi Yin (ziyi.yin@gatech.edu)
## ruined by: Mark Glines (mark@extreme-scale.com)

using LoggingExtras, Dates

const date_format = "yyyy-mm-dd HH:MM:SS"

timestamp_logger(logger) = TransformerLogger(logger) do log
  merge(log, (; message = "$(Dates.format(now(), date_format)) $(log.message)"))
end

ConsoleLogger(stdout, Logging.Info) |> timestamp_logger |> global_logger

@info("loading libraries")

using DrWatson
@quickactivate "proxy-example"
using Pkg; Pkg.instantiate()
using Plots
using SegyIO
using Statistics
using Polynomials
using JutulDarcyRules
include(srcdir("utils.jl"))

if ~isdir("slim.gatech.edu/data//synth/Compass")
    @info("pre-fetching data files")
    run(`wget -r ftp://slim.gatech.edu/data//synth/Compass`) # this might take a while
end

@info("reading velocity tensor")

# get velocity
block = segy_read("slim.gatech.edu/data/synth/Compass/final_velocity_model_ieee_6m.sgy")

# original compass model is a grid of 1911x2730x341 cells of size 25m*25m*6m
n_orig = (1911,2730,341)
d_orig = (25f0,25f0,6f0)

@info("reading cell x,y positions")

sx = get_header(block, "SourceX")
sy = get_header(block, "SourceY")

@info("converting data to 3d array")

# turn it to 3D cube
v = zeros(Float32,n_orig)
for i = 1:n_orig[1]
    x = d_orig[1].*(i-1)
    inds = findall(sx.==x)
    slice = block.data[:,inds[sortperm(sy[inds])]]

    v[i,:,:] = transpose(slice)/1f3
end

seed = 2
rng = Random.seed!(seed)

@info("extracting chunk at random position using seed $seed")

# take a 4km chunk of input
top_layer = 182
chunk_shape_xy_km = (4, 4)
chunk_shape = (floor(Int, chunk_shape_xy_km[1]*1000/d_orig[1]), floor(Int, chunk_shape_xy_km[2]*1000/d_orig[2]), n_orig[end] - top_layer + 1)
xoff = 1
yoff = 1 + floor(Int, n_orig[2]/2)
println("xoff is $xoff, yoff is $yoff")
chunk = v[xoff:xoff+chunk_shape[1]-1, yoff:yoff+chunk_shape[2]-1, top_layer:end]

@info("downsampling velocity array from $(size(chunk))")

factor = (1,1,1)                                                            # the full 160³
# factor = (2,2,2)                                                            # downsampling to 80³
# factor = (8,8,8)                                                            # downsampling to 20³
if factor == (1,1,1)
    vsmall = Float64.(chunk)
else
    vsmall = Float64.(1f0./downsample(1f0./chunk, factor))
end
h = (top_layer-1) * 1.0 * d_orig[end]                                            # starting depth of the simulation model

@info("converting $(size(vsmall)) velocity array to permeabilities")

Kh = VtoK_stochastic.(rng, maximum(vsmall), vsmall) # horizontal permeability
# Kh = VtoK_simple.(vsmall) # horizontal permeability
# Kh = VtoK_banded.(vsmall) # horizontal permeability
K = Kh * md                 # millidarcy
ϕ = Ktoϕ.(Kh)               # porosity

kvoverkh = 0.1                 # kv/kh ratio
n = size(K)                    # model size for flow simulation
d = Float64.(d_orig .* factor) # discretization

xmid = Int(floor(n[1] * 0.5))
ymid = Int(floor(n[2] * 0.5))
zvec = zvec = Kh[xmid, ymid, :]                 # find the point in the center column with the highest permeability
zstart = Int(floor(n[3] * 0.65))                # between 65% of the way down
zend   = Int(floor(n[3] * 0.85))                # and 85% of the way down
zstart = argmax(zvec[zstart:zend]) + zstart - 2 # and use that as our injection point
zend = zstart + 2

tstep = 365.25 * 2 * ones(30)     # 60 years, in 2-year increments
tstepidle = 365.25 * 2 * ones(25) # 50 years, in 2-year increments

@info("setting up jutul simulation")

model = jutulModel(n, d, vec(padϕ(ϕ)), K1to3(K; kvoverkh=kvoverkh), h)

## simulation time steppings
tot_time = sum(tstep)

## injection well location
pore_volumes = sum(ϕ[2:end-1,2:end-1,1:end-1] .* (vsmall[2:end-1,2:end-1,1:end-1].>3.5)) * prod(d) ## reservoir pore volumes (v>3.5 is reservoir)
storage_capacity = 0.05                                                                            ## 20% for 2D, 5% for 3D
irate = storage_capacity * pore_volumes / tot_time / 24 / 60 / 60                                  ## set injection rate for storage capacity
q   = jutulVWell(irate, (xmid * d[1], ymid * d[2]); startz = zstart * d[end], endz = zend * d[end])  ## well with rates and locations
qidle = jutulVWell(0.0, (xmid * d[1], ymid * d[2]); startz = zstart * d[end], endz = zend * d[end])  ## inactive well
## set up modeling operator
S     = jutulModeling(model, tstep)
Sidle = jutulModeling(model, tstepidle)

## simulation
mesh_ = CartesianMesh(model)            ## cartesian mesh
T(x) = log.(KtoTrans(mesh_, K1to3(exp.(x); kvoverkh=kvoverkh)))         ## convert log permeability to log transmissibility (which jutul prefers)

## log permeability
logK = log.(K)

@info("running jutul simulation")

## simulation
@time state     = S(    T(logK), vec(padϕ(ϕ)), q)
@time stateidle = Sidle(T(logK), vec(padϕ(ϕ)), qidle; state0=state.states[end])

@info("plotting results")

wellx = []
welly = []
wellz = []
for z in zstart:zend
    push!(wellx, xmid)
    push!(welly, ymid)
    push!(wellz, z)
end

daylen = Int(floor(sum(tstep)+sum(tstepidle)))
daylen = sizeof("$daylen")
days = 0

states = vcat(state.states, stateidle.states)
tstepall = vcat(tstep, tstepidle)
for (timestep, thisstate) in zip(tstepall, states)
    global days
    days += timestep
    sat = reshape(Saturations(thisstate), n...)
    pre = reshape(Pressure(thisstate   ), n...)
    perm = exp.(logK)/md
    maxperm = maximum(perm)
    xsize = 2
    ### plotting xz planes
    xzvel = heatmap(perm[:,ymid,:]', seriescolor=:ice, yflip=true, clims=(100, maxperm), colorbar_scale=:log10, colorbar_title="m²", title="true permeability (XZ)")
    scatter!(wellx, wellz, markershape=:x, markercolor=:yellow, markersize=xsize, legend=false)
    xzsat = heatmap(sat[:,ymid,:]', seriescolor=:ice, yflip=true, colorbar_title="% pore volume", title="saturation (XZ)")
    scatter!(wellx, wellz, markershape=:x, markercolor=:black, markersize=xsize, legend=false)
    xzpre = heatmap((pre[:,ymid,:]') ./ 1000000, seriescolor=:ice, yflip=true, colorbar_title="10⁶ Pascal", title="pressure (XZ)")
    scatter!(wellx, wellz, markershape=:x, markercolor=:black, markersize=xsize, legend=false)
    ### plotting yz planes
    yzvel = heatmap(perm[xmid,:,:]', seriescolor=:ice, yflip=true, clims=(100, maxperm), colorbar_scale=:log10, colorbar_title="m²", title="true permeability (YZ)")
    scatter!(wellx, wellz, markershape=:x, markercolor=:yellow, markersize=xsize, legend=false)
    yzsat = heatmap(sat[xmid,:,:]', seriescolor=:ice, yflip=true, colorbar_title="% pore volume", title="saturation (YZ)")
    scatter!(wellx, wellz, markershape=:x, markercolor=:black, markersize=xsize, legend=false)
    yzpre = heatmap((pre[xmid,:,:]') ./ 1000000, seriescolor=:ice, yflip=true, colorbar_title="10⁶ Pascal", title="pressure (YZ)")
    scatter!(wellx, wellz, markershape=:x, markercolor=:black, markersize=xsize, legend=false)
    l = @layout [a b; c d; e f]
    plot(xzvel, yzvel, xzsat, yzsat, xzpre, yzpre, dpi=200, size=(1080,720), layout=l)
    daystring = lpad(Int(floor(days)), daylen, "0")
    savefig(plotsdir("after$(daystring)days.png"))
end
