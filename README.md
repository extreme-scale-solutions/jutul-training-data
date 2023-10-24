# jutul-training-data

This project runs some CO₂ injection simulations, using the Julia Jutul library.
It is heavily based on [Zizi Yin's proxy example code](https://github.com/slimgroup/proxy-example).

# Pieces

There are a few moving pieces to this, which live in different places.

## COMPASS data model

The slimgroup folks have shared their compass dataset online, and we use random subsets of it for our simulations.

If the data is not present, it will be downloaded.  It's pretty big, so it's easier to just use the copy we already have.

Create a symlink in your jutul-training-data folder, to point to the copy in CFS:

```
ln -s /global/cfs/cdirs/m3863/mark/slim.gatech.edu .
```

## Physical system

All simulation data is based on subsets of the Compass model, which is a grid of 1911x2730x341 cells of size 25m\*25m\*6m.
We select a random 2 km² patch of land from this model, and take that subset with a fixed depth, producing a 160³ tensor.
Each cell in this tensor is a chunk of land, whose size is 25\*25\*6 meters.
This 160³ tensor can then be simulated as-is (which takes about 1 hour), or downsampled to 80³, 40³ or 20³ and then we
simulate the downsampled system.  The smaller scales are very useful for development and testing.

Each actual physics simulation lives in its own data folder.  A data folder might have a path like `data/v5/20³/1016296115306602175`.
* "v5" is the version of the data format; it's evolved quite a bit to make sure all of the necessary fields are present.
* "20³" is the input tensor shape: (20,20,20).
* "1016296115306602175" is the random seed used to choose this simulation's position (within the Compass dataset), and the random variations on the input permeability tensor

Within this folder are the following files:
* inputs.jld2 contains all of the input parameters and tensors.
* meta.json contains the same data as inputs.jld2, minus the big tensors.

The script `mkinput.jl` generates one or more of these, for a given downsampling rate.  The parameters are the parent folder to create seed folders underneath, the downsampling rate, and the number of seeds to create.

```
julia --project scripts/mkinput.jl ~/PSCRATCH/v5/20³  8 4000
julia --project scripts/mkinput.jl ~/PSCRATCH/v5/40³  4 4000
julia --project scripts/mkinput.jl ~/PSCRATCH/v5/80³  2 4000
julia --project scripts/mkinput.jl ~/PSCRATCH/v5/160³ 1 4000
```

This generates 4000 inputs for each data size, and puts each in their own dataset folder.

## Physics simulation

The script `simulate-batch.jl` runs the actual physics simulation on one or more input folders, and adds these files to the existing data dir:

* outputs.jld2 contains all of the stuff in meta.json, plus arrays of saturation and pressure output tensors.
* after\*days.png: renderings of the system state at various timestamps
* $seed.mp4: a video file of the after\*days.png frames

```
sbatch --account=$ACCOUNT --constraint=cpu --nodes=1 --qos=regular --time=12:00:00 ./simulate-batch-10.sh ~/PSCRATCH/v5/20³
sbatch --account=$ACCOUNT --constraint=cpu --nodes=1 --qos=regular --time=12:00:00 ./simulate-batch-10.sh ~/PSCRATCH/v5/40³
sbatch --account=$ACCOUNT --constraint=cpu --nodes=1 --qos=regular --time=12:00:00 ./simulate-batch-10.sh ~/PSCRATCH/v5/80³
sbatch --account=$ACCOUNT --constraint=cpu --nodes=1 --qos=regular --time=12:00:00 ./simulate-batch-6.sh ~/PSCRATCH/v5/160³
```

Make sure ACCOUNT is set to your billing account id; ours is m3863.

We run 10 simulations in parallel for most of the data sizes, as this is a good level of parallelism for this system.
We run 6 simulations of the largest scale, so that they won't run out of memory.


# Other commands

These are commands I found useful while generating the full dataset.

```
# set up julia environment
module use $CFS/m3863/mark/lmod
module load julia/1.8.5-mark
# see how many finished simulations are in cfs
cd ~/CFS/training-data/training-samples/v5; for SIZE in 20 40 80 160; do ls ${SIZE}³ | wc -l; done
# spawn a simulation job for a specified size
cd ~/src/proxy-example; sbatch --account=m3863 --constraint=cpu --nodes=1 --qos=regular --time=12:00:00 ./simulate-batch-10.sh ~/PSCRATCH/v5/20³
cd ~/src/proxy-example; sbatch --account=m3863 --constraint=cpu --nodes=1 --qos=regular --time=12:00:00 ./simulate-batch-10.sh ~/PSCRATCH/v5/40³
cd ~/src/proxy-example; sbatch --account=m3863 --constraint=cpu --nodes=1 --qos=regular --time=12:00:00 ./simulate-batch-10.sh ~/PSCRATCH/v5/80³
cd ~/src/proxy-example; sbatch --account=m3863 --constraint=cpu --nodes=1 --qos=regular --time=12:00:00 ./simulate-batch-6.sh ~/PSCRATCH/v5/160³
# move finished things from pscratch to cfs
cd ~/PSCRATCH/v5; for SIZE in 20 40 80 160; do echo SIZE=$SIZE; cd ${SIZE}³; ls | while read THING; do if test -f $THING/after40177days.png; then echo $THING; mv $THING /global/homes/i/infinoid/CFS/training-data/training-samples/v5/${SIZE}³; fi; done; cd ..; done
```

# Troubleshooting

## CUDA

These commands may generate an error message with a backtrace, from the CUDA library not finding any GPUs (on CPU nodes).

```
Error: 2023-05-31 06:17:11 Failed to initialize CUDA
```

This message and backtrace can be safely ignored.

## GKS/Qt/xcb

There's an `export GKSwstype=nul` in the simulate-batch shell scripts.  This reduces confusion in Plots/GR libraries when we run in headless (no GUI) mode.

If you see a bunch of errors like this, during the plot phase:

```
connect: Connection refused
GKS: can't connect to GKS socket application

GKS: Open failed in routine OPEN_WS
GKS: GKS not in proper state. GKS must be either in the state WSOP or WSAC in routine ACTIVATE_WS
qt.qpa.xcb: could not connect to display
qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.

Available platform plugins are: linuxfb, minimal, offscreen, vnc, xcb.
```

It's because you didn't set this variable.  This is not fatal… it's just spammy and slow.
