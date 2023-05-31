### Read the Compass slice from SLIM website and take out the 2D slice
## Please refer to the details in the MIT license of this repository and in the license of the Compass model
## author: Ziyi Yin (ziyi.yin@gatech.edu)
## ruined by: Mark Glines (mark@extreme-scale.com)

using LoggingExtras, Dates

if length(ARGS) < 3
    println("Usage: julia --project scripts/mkinput.jl <datasetdir> <divisor> <num_inputs>")
    exit(1)
end
datasetdir = ARGS[1]
divisor = parse(Int, ARGS[2])
num_inputs = parse(Int, ARGS[3])

const date_format = "yyyy-mm-dd HH:MM:SS"

timestamp_logger(logger) = TransformerLogger(logger) do log
  merge(log, (; message = "$(Dates.format(now(), date_format)) $(log.message)"))
end

ConsoleLogger(stdout, Logging.Info) |> timestamp_logger |> global_logger

@info("loading libraries")

using DrWatson
@quickactivate "proxy-example"
using Pkg; Pkg.instantiate()
using SegyIO
using Statistics
using Polynomials
using JutulDarcyRules
import Random
import JLD2
import JSON

include(srcdir("utils.jl"))

@info("caching data files")

if ~isdir("slim.gatech.edu/data//synth/Compass")
    run(`wget -r ftp://slim.gatech.edu/data//synth/Compass`) # this might take a while
end

@info("reading Compass velocity model")

# get velocity
block = segy_read("slim.gatech.edu/data/synth/Compass/final_velocity_model_ieee_6m.sgy")
# block is a matrix of size 341x5217030
# that is Z x (X*Y)

# original compass model is a grid of 1911x2730x341 cells of size 25m*25m*6m
n_orig = (1911,2730,341)
d_orig = (25f0,25f0,6f0)

@info("reading grid coordinate tables")

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

for iter in 1:num_inputs
    @info("generating random seed")

    seed = abs(rand(Int))

    sampledir = "$datasetdir/$seed"

    try
        mkdir(sampledir)
    catch
        @warn "skipping duplicate seed $seed"
        continue
    end

    rng = Random.seed!(seed)

    @info("extracting chunk at random position using seed $seed")

    # take a 4km chunk of input
    top_layer = 182
    chunk_shape_xy_km = (4, 4)
    chunk_shape = (floor(Int, chunk_shape_xy_km[1]*1000/d_orig[1]), floor(Int, chunk_shape_xy_km[2]*1000/d_orig[2]), n_orig[end] - top_layer + 1)
    xoff = abs(rand(rng, Int)) % (n_orig[1]-chunk_shape[1]-1) + 1
    yoff = abs(rand(rng, Int)) % (n_orig[2]-chunk_shape[2]-1) + 1
    println("xoff is $xoff, yoff is $yoff")
    chunk = v[xoff:xoff+chunk_shape[1]-1, yoff:yoff+chunk_shape[2]-1, top_layer:end]

    @info("downsampling velocity array from $(size(chunk))")

    # factor = (1,1,1)                                                            # the full 160³
    # factor = (2,2,2)                                                            # downsampling to 80³
    # factor = (4,4,4)                                                            # downsampling to 40³
    # factor = (8,8,8)                                                            # downsampling to 20³
    factor = (divisor, divisor, divisor)
    if factor == (1,1,1)
        vsmall = Float64.(chunk)
    else
        vsmall = Float64.(1f0./downsample(1f0./chunk, factor))
    end
    h = (top_layer-1) * 1.0 * d_orig[end]                                            # starting depth of the simulation model

    @info("converting $(size(vsmall)) velocity array to permeabilities")

    Kh = VtoK_stochastic.(rng, maximum(vsmall), vsmall)          # horizontal permeability
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

    tot_time = sum(tstep)
    pore_volumes = sum(ϕ[2:end-1,2:end-1,1:end-1] .* (vsmall[2:end-1,2:end-1,1:end-1].>3.5)) * prod(d)      ## reservoir pore volumes (v>3.5 is reservoir)
    storage_capacity = 0.05                                                                                 ## 20% for 2D, 5% for 3D
    irate = storage_capacity * pore_volumes / tot_time / 24 / 60 / 60                                       ## set injection rate for storage capacity

    @info("saving input parameters")

    savedict = Dict(
        "n"         => n,
        "d"         => d,
        "kvoverkh"  => kvoverkh,
        "factor"    => factor,
        "seed"      => seed,
        "h"         => h,
        "xmid"      => xmid,
        "ymid"      => ymid,
        "xoff"      => xoff,
        "yoff"      => yoff,
        "irate"     => irate,
        "zstart"    => zstart,
        "zend"      => zend,
        "tstep"     => tstep,
        "tstepidle" => tstepidle,
    )
    open("$sampledir/meta.json", "w") do f
        redirect_stdout(f) do
            JSON.print(savedict, 2)
        end
    end
    savedict["K"] = K
    savedict["ϕ"] = ϕ
    savedict["vsmall"] = vsmall
    JLD2.save("$sampledir/inputs.jld2", savedict)
end
