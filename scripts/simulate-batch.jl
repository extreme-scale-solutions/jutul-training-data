### Read the Compass slice from SLIM website and take out the 2D slice
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
using Statistics
using JutulDarcyRules
import Random
import JLD2
import JSON

include(srcdir("utils.jl"))

@info("ARGS is $ARGS")
for sampledir in ARGS
    @info("reading inputs from $sampledir")

    savedict  = JLD2.load("$sampledir/inputs.jld2")
    d         = savedict["d"]
    ϕ         = savedict["ϕ"]
    kvoverkh  = savedict["kvoverkh"]
    factor    = savedict["factor"]
    K         = savedict["K"]
    n         = savedict["n"]
    h         = savedict["h"]
    vsmall    = savedict["vsmall"]
    xmid      = savedict["xmid"]
    ymid      = savedict["ymid"]
    xoff      = savedict["xoff"]
    yoff      = savedict["yoff"]
    zstart    = savedict["zstart"]
    zend      = savedict["zend"]
    tstep     = savedict["tstep"]
    tstepidle = savedict["tstepidle"]

    @info("setting up jutul simulation")

    model = jutulModel(n, d, vec(padϕ(ϕ)), K1to3(K; kvoverkh=kvoverkh), h)

    ## simulation time steppings
    tot_time = sum(tstep)

    ## injection well location
    pore_volumes = sum(ϕ[2:end-1,2:end-1,1:end-1] .* (vsmall[2:end-1,2:end-1,1:end-1].>3.5)) * prod(d)      ## reservoir pore volumes (v>3.5 is reservoir)
    storage_capacity = 0.05                                                                                 ## 20% for 2D, 5% for 3D
    irate = storage_capacity * pore_volumes / tot_time / 24 / 60 / 60                                       ## set injection rate for storage capacity
    q   = jutulVWell(irate, (xmid * d[1], ymid * d[2]); startz = zstart * d[end], endz = zend * d[end])       ## well with rates and locations
    qidle = jutulVWell(0.0, (xmid * d[1], ymid * d[2]); startz = zstart * d[end], endz = zend * d[end])       ## inactive well

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
    @time state         = S(T(logK), vec(padϕ(ϕ)), q)
    @time stateidle = Sidle(T(logK), vec(padϕ(ϕ)), qidle; state0=state.states[end])

    @info("saving outputs")
    states = vcat(state.states, stateidle.states)
    savedict = Dict(
        "d"           => d,
        "kvoverkh"    => kvoverkh,
        "factor"      => factor,
        "n"           => n,
        "h"           => h,
        "xoff"        => xoff,
        "yoff"        => yoff,
        "xmid"        => xmid,
        "ymid"        => ymid,
        "zstart"      => zstart,
        "zend"        => zend,
        "tstep"       => tstep,
        "irate"       => irate,
        "tstepidle"   => tstepidle,
        "numstates"   => size(state.states)[1],
        "saturations" => map(x -> reshape(Saturations(x), n...), states),
        "pressures"   => map(x -> reshape(Pressure(   x), n...), states),
    )
    JLD2.save("$sampledir/outputs.jld2", savedict)

    @info("plotting results")

    wellx = []
    welly = []
    wellz = []
    for z in zstart:zend
        push!(wellx, xmid)
        push!(welly, ymid)
        push!(wellz, z)
    end

    daylen = Int(floor(sum(tstep) + sum(tstepidle)))
    daylen = sizeof("$daylen")
    days = 0
    tstepall = vcat(tstep, tstepidle)
    for (timestep, thisstate) in zip(tstepall, states)
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
        savefig("$sampledir/after$(daystring)days.png")
    end
end
