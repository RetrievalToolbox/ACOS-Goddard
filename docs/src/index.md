# ACOS Goddard

Welcome to ACOS Goddard, an independent implementation of the ACOS algorithm that was originally created at NASA JPL. The NASA JPL implementation of the algorithm is `RtRetrievalFramework` and is available [here](https://github.com/nasa/RtRetrievalFramework). Our implementation of ACOS uses the `RetrievalToolbox` library, which can be freely downloaded [here](https://www.github.com/US-GHG-Center/RetrievalToolbox.jl). One can think of this as an own application that is built on top of `RetrievalToolbox`.

The purpose of this application is to serve as a demonstration for science users to see how a moderately complex CO$_2$ retrieval algorithm can be built and run with the `RetrievalToolbox` software library.

Before downloading or using this application, users are highly encouraged to read the [important notes](important.md). Further, there are notable differences to the reference implementation as well as known issues, which are listed [here](known_issues.md).

## Installation

### Installing Julia and RetrievalToolbox

We recommend installing Julia via JuliaUp. Instructions can be found [here](https://github.com/JuliaLang/juliaup) for various platforms. RetrievalToolbox can be directly installed via Julia's package manager, however the XRTM radiative transfer solver library must be installed manually. Instructions for both are found in the RetrievalToolbox [repository](https://www.github.com/US-GHG-Center/RetrievalToolbox.jl]).

### Installing packages needed for ACOS-Goddard

The most convenient way to get started is to jump into the directory in which you have cloned ACOS-Goddard and install all required packages automatically. This is done by invoking Julia with the `project` argument:

``` bash
# (within the cloned directory)
julia --project=./
```

And from here, let Julia install all required packages:
``` julia
using Pkg;
Pkg.add(path="https://www.github.com/US-GHG-Center/RetrievalToolbox.jl");
Pkg.resolve();
Pkg.instantiate();
```

This should ideally install the needed modules to run the example provided here.

### Obtaining required auxiliary data

In order to run the provided example(s), users must download additional data not present in this repository. First, a solar model is required. Our ACOS implementation currently supports two different solar model formats. One is the solar model file that is used in RtRetrievalFramework, which is the best choice for users who want to stick closer to the operational OCO code. That solar model can be downloaded from the Github repository with the following shell commands

``` bash
cd example_data
wget https://github.com/nasa/RtRetrievalFramework/raw/refs/heads/master/input/common/input/l2_solar_model.h5
```

Alternatively, users can make use of the TSIS-1 Hybrid Solar Reference Spectrum, developed and published by the Laboratory for Atmospheric and Space Physics (LASP) at the University of Colorado Boulder. The download portal is found [here](https://lasp.colorado.edu/lisird/data/tsis1_hsrs_p1nm), or again use the following shell command to directly download the model file at 0.001 nm spectral resolution:

``` bash
# (assuming you are still in `example_data`)
wget https://lasp.colorado.edu/lisird/resources/lasp/hsrs/v2/hybrid_reference_spectrum_p005nm_resolution_c2022-11-30_with_unc.nc
```

Note that an appropriate change in the corresponding command line argument to the application must be done. The "--solar_model" command line argument must point to the right file. No further changes are needed, the application detects the type of model based on the file ending (`.nc` for the LASP TSIS file, and `.h5` for the RtRetrievalFramework format). By default, examples use the model as provided by RtRetrievalFramework.

Further, the aerosol optical properties as well as some other additional data (such as the non-uniform sampling grid) are required:

``` bash
# (assuming you are still in `example_data`)
wget https://github.com/nasa/RtRetrievalFramework/raw/refs/heads/master/input/common/input/l2_aerosol_combined.h5
wget https://github.com/nasa/RtRetrievalFramework/raw/refs/heads/master/input/oco/input/l2_oco_static_input.h5
```

Note that the location of these two files is *hard-coded* inside our application (`./example_data/...`), so if you want to replace these files or their location, you have to manually make the appropriate changes inside the `*.jl` files.

### Obtaining Spectroscopy

As of now, RetrievalToolbox supports spectroscopy in the ABSCO format only. To use this retrieval algorithm, users need to obtain the spectroscopy tables for O$_2$, CO$_2$ and H$_2$O. The GES-DISC hosts a [webpage](https://disc.gsfc.nasa.gov/information/glossary?title=OCO-2%20ABSCO) for the time being which lists contact addresses where researchers can request those files.


## Running Examples

Running the provided example retrieval is done as follows, assuming that `XRTM_PATH` is correctly set and the needed auxiliary files are downloaded and placed in the correct directories. It is important that the `--project=./` argument is supplied, otherwise the Julia session has no knowledge of the additional packages needed by this application!

``` bash
XRTM_PROGRESS=1 julia --project=./ ./from_repl.jl
```

or from within Julia, when started with `--project=./`
``` julia
include("from_repl.jl")
```

The `from_repl.jl` script includes some reasonable default settings to run the provided example retrieval.

Note that RetrievalToolbox supports multi-threading for the monochromatic RT calculations, and the number of threads that can be utilized is determined by the `JULIA_NUM_THREADS` environment variable. So to speed up the demo retrieval, use

``` bash
XRTM_PROGRESS=1 JULIA_NUM_THREADS=3 julia --project=./ ./from_repl.jl
```

For this particular set-up, choosing `3` threads tends to provide the best benefit, using more threads does not achieve notable speed-up in our experience. The high-accuracy approximation technique (LSI) cannot yet benefit from multi-threading. For parallel processing of multiple **soundings** via, e.g., GNU Parallel, users might want to benchmark which combination of parallel instances vs. the number of threads to use. The memory footprint also changes when using more than one thread.

For new users, we recommend starting with the interactive notebooks inside the `examples` folder. A working Jupyter installation is required, and it must be able to access the Julia kernel(s). In order to register Julia kernels with Jupyter (or JupyterLab), simply install `IJulia` by typing (MacOS/Linux) `julia -e 'using Pkg; Pkg.add("IJulia")'` or manually inside the Julia REPL:

``` julia
] add IJulia
```

Once `IJulia` is installed through the package manager, a local notebook session can be spawned by typing

``` julia
notebook()
```

which will start the Jupyter server and open a web browser. Alternatively, users can also spawn Jupyter(Lab) from their preferred shell (`jupyter`, or `juptyer-lab`) and then navigate to the location of the example notebooks.

### Command line arguments for the `run.jl` script

Instead of running `from_repl.jl`, users can of course tailor the inputs to the retrieval algorithm. Below is an example of how to call `run.jl`, which is the main application to be used. Like before, this assumes that the required auxiliary data is downloaded and placed into the correct folders, and that `XRTM_PATH` is correctly set. The below example runs a single retrieval with the sounding ID `2021030111564431`. To process several soundings, users can invoke the `--sounding_id_list` command line argument, explained further below.

``` bash
XRTM_PROGRESS=1 JULIA_NUM_THREADS=3 julia --project=./ ./run.jl \
        --solar_model ./example_data/l2_solar_model.h5 \
        --L1b ./example_data/2021030111564431_inputs.h5 \
        --L2Met ./example_data/2021030111564431_inputs.h5 \
        --L2CPr ./example_data/2021030111564431_inputs.h5 \
        --output 2021030111564431.h5 \
        --o2_spec ./example_data/o2_v52.hdf \
        --o2_scale 1.0048 \
        --co2_spec ./example_data/co2_v52.hdf \
        --co2_scale_weak 0.994 \
        --co2_scale_strong 0.998 \
        --h2o_spec ./example_data/h2o_v52.hdf \
        --sounding_id 2021030111564431 \ # alternative: --sounding_id_list sounding_id_list.txt
        --spec 1,2,3 \
        --polarized true \
        --Lambertian false \
        --aerosols true \
        --retrieve_aerosols true \
        --LSI true \
        --Nhigh 16 \
        --dsigma_scale 2.0 \
        --gamma 50.0 \
        --max_iterations 10
```

Each command line argument will be explained below.

* `--solar_model`: Points to the path of the solar model file, either the one that ships with RtRetrievalFramework, or the LASP TSIS-1 hybrid reference spectrum (see above on how to download them).

* `--L1b`: Points to the path of an OCO-2/3 L1b file (calibrated, geo-located measurements). Must be downloaded from the GES-DISC: [OCO-2/3 L1B collection](https://disc.gsfc.nasa.gov/datasets?keywords=oco2%20oco3&page=1&processingLevel=1B)

* `--L2Met`: Points to the path of an OCO-2/3 L2Met file (meteorology and aerosols sampled at measurement locations). Must be downloaded from GES-DISC: [OCO-2 L2Met collection](https://disc.gsfc.nasa.gov/datasets?keywords=OCO2_L2_Met_11.2r&page=1) or [OCO-3 L2Met collection](https://disc.gsfc.nasa.gov/datasets?keywords=OCO3_L2_Met_11r&page=1)

* `--L2CPr`: Points to the path of an OCO-2/3 L2CPr file (CO2 prior profiles at measurement locations). Must be downloaded from the GES-DISC: [OCO-2 L2CPr collection](https://disc.gsfc.nasa.gov/datasets?keywords=OCO2_L2_CO2Prior_11.2r&page=1) or [OCO-3 L2CPr collection](https://disc.gsfc.nasa.gov/datasets?keywords=OCO3_L2_CO2Prior_11r&page=1) NOTE! The three files (L1b, L2Met and L2CPr) must belong to the same orbit to contain the data associated with one single **sounding ID**!

* `--sounding_id`: Must be a typical 16-digit OCO-type sounding ID, which identifies the single sounding to be processed.

* `--sounding_id_list`: Points to the path of a file that contains a list of 16-digit OCO-type sounding IDs. The file should be a text file with sounding IDs, one per line. When performing retrievals for many scenes, it can be advantageous to use this open instead of calling the algorithm for every scene individually, as one suffers the compilation step every time the algorithm is freshly called.

* `--output`: Path to the directory in which the output file will be written. **NOTE!** Within this application, any existing output file is automatically **overwritten**. If you need to skip existing files, it is best to implement those checks within the script(s) that execute the application. Alternatively, users can modify the code in `main.jl`. A single output file is created for each individual sounding ID; users who need a different solution (e.g. writing results into a big parallel NetCDF file) will need to write their own code.

* `--touch_init`: (optional) Determines whether the program will write an *init*-file to signify that this particular sounding ID has started processing. This is useful for book-keeping and to isolate soundings that make the program crash. Once the output file is successfully written out, the *init*-file is removed. The default is `true`.

* `--o2_spec`: Points to the path of the oxygen spectroscopy ABSCO file.

* `--o2_scale`: (optional) Must be a real number. The entire oxygen spectroscopy table is scaled by this value. For ABSCO version 5.2, the recommended value is `1.0048`, but that number was derived with the RtRetrievalFramework and may not be the ideal choice for this algorithm. The default is `1.0`.

* `--co2_spec`: Points to the path of the carbon dioxide spectroscopy ABSCO file.

* `--co2_scale_weak`: (optional). Must be a real number. The "weak" part of the CO2 spectroscopy table (< 2 µm) is is scaled by this number. For ABSCO version 5.2, the recommended value is `0.994`, but that number was derived with the RtRetrievalFramework and may not be the ideal choice for this algorithm. The default is `1.0`.

* `--co2_scale_strong`: (optional). Must be a real number. The "strong" part of the CO2 spectroscopy table (> 2 µm) is is scaled by this number. For ABSCO version 5.2, the recommended value is `0.998`, but that number was derived with the RtRetrievalFramework and may not be the ideal choice for this algorithm. The default is `1.0`.

* `--h2o_spec`: Points to the path of the water vapor spectroscopy ABSCO file.

* `--spec`: A comma-separated list of integers (1,2 and/or 3) that represent which spectrometers (1 = O$_2$ A-band, 2 = CO$_2$ band at 1.6 µm, 3 = CO$_2$ band at 2.06 µm) are used in the retrieval. This can be something like `1,3`, `2` (for a single-band retrieval of CO$_2$), or to use all (standard) `1,2,3`. The spectrometer also impacts the state vector - only those state vector elements will be used which make sense given the spectrometer configuration. For example, the CO$_2$ profile will only be retrieved if either spectrometer `2` or `3` are used. Surface pressure and SIF are not retrieved if spectrometer `1` is not in the list.

* `--polarized`: (optional) Must be `true` or `false`. If set to `true`, polarization is accounted for in the model, and the appropriate settings are chosen for the radiative transfer. The algorithm is faster if this is set to `false`, however note that accounting for polarization is required for accurate CO$_2$ retrievals from OCO-2/3 due to the polarization sensitivity of the instrument. The default is `true`.

* `--Lambertian`: (optional) Must be `true` or `false`. If set to `true`, a Lambertian BRDF is used as the surface model in the radiative transfer calculations. If set to `false`, the Rahman-Pinty-Verstraete (RPV) type BRDF is used, as laid out in the ACOS ATBD. The default is `false`.

* `--aerosols`: (optional) Must be `true` or `false`. If set to `true`, aerosols are part of the model atmosphere. The types and abundances are selected via the L2Met file. As in ACOS, the two most abundant tropospheric types are chosen, along with a high-altitude stratospheric sulphur type aerosol, a water cloud type and an ice cloud type. If set to `false`, none of these aerosols are considered in the algorithm at all, which results in significant speed-up of the run time (~5s per iteration). Aerosols are generally required for accurate retrievals. The default is `true`.

* `--retrieve_aerosols`: (optional) Must be `true` or `false`. If set to `true`, all aerosols that are in the atmosphere are being retrieved with all three shape parameters (total AOD, height and width of the aerosol layer). Setting this to `false` will still retain the aerosol layers in the model atmosphere, however none are being adjusted during the retrieval and will keep their first-guess value as informed by the L2Met (which samples the GEOS-IT model). The default is `true`.

* `--LSI`: (optional) Must be `true` or `false`. If set to `true`, the algorithm will use the *Low Streams Interpolation* technique to approximate high-accuracy multiple scattering RT calculations. If set to `false`, no high-accuracy RT computations are made and only the single-scattering and two-stream RT computations are performed. The default is `true`.

* `--Nhigh`: (optional) Must be an even integer >= 2. If `LSI` is set, `Nhigh` will the number of (full) quadrature streams used for the high-accuracy approximation. Larger numbers mean higher accuracy, but slower runtime. The default value is `16`.

* `--dsigma_scale`: (optional) Must be a real number. This is the scale factor applied to the $d\sigma^2$ calculation used to establish convergence of the retrieval (see OCO-2 ATBD Section 3.5.1). The default value is `2.0`.

* `--gamma`: (optional) Must be a real number. This is the Levenberg-Marquardt $\gamma$ parameter that controls the step size of the inversion. The default value is `10.0`, however we recommend higher values (e.g. `500`) as our algorithm is known to behave more non-linearly when getting closer to the solution.

* `--max_iterations`: (optional) Must be an integer. The inversion is halted after `--max_iterations` number of iterations. Setting this to `1`, for example, will cause the algorithm to execute the forward model once, update the state vector once, and then quit. Setting this to `0` will skip the forward model execution, but still return the objects that could be used for introspection of the retrieval scene setup when used interactively. The default value is `10`.

## Parallel processing

The application supports rudimentary parallel processing via Julia's built-in functions for distributed computing. From a user perspective, the only required change is to execute the script with multiple processes (or workers) with the `-p` flag, and ACOS Goddard will take care of the rest. In our example shipped in this repository, one would run 3 parallel processes in total to process all 3 IDs in the `sounding_id_list.txt` files at the same time. Note that the `-p` flag denotes the number of additional worker processes, so the total number of available processes is that number plus one.

``` bash
XRTM_PROGRESS=1 JULIA_NUM_THREADS=1 julia --project=./ -p 2 ./run.jl \
        --solar_model ./example_data/l2_solar_model.h5 \
        --L1b ./example_data/2021030111564431_inputs.h5 \
        --L2Met ./example_data/2021030111564431_inputs.h5 \
        --L2CPr ./example_data/2021030111564431_inputs.h5 \
        --output 2021030111564431.h5 \
        --touch_init true \
        --o2_spec ./example_data/o2_v52.hdf \
        --o2_scale 1.0048 \
        --co2_spec ./example_data/co2_v52.hdf \
        --co2_scale_weak 0.994 \
        --co2_scale_strong 0.998 \
        --h2o_spec ./example_data/h2o_v52.hdf \
        --sounding_id_list sounding_id_list.txt \
        --spec 1,2,3 \
        --polarized true \
        --Lambertian false \
        --aerosols true \
        --retrieve_aerosols true \
        --LSI true \
        --Nhigh 16 \
        --dsigma_scale 2.0 \
        --gamma 50.0 \
        --max_iterations 10
```

Note that the ACOS Goddard application is not only spawned three times - we are making use of the `SharedArray` type for distributed computing in Julia such that the big spectroscopy tables are only loaded into memory **once**, and all worker processes access the same data. The `-p` flag determines the number of **additional** processes that are spawned; so `-p 2` will spawn 2 more processes, and ACOS Goddard will perform retrievals on a total of 3 parallel processes. In multi-processing mode, logging outputs are filtered, so any logging message that does not start with `[MAIN]` is suppressed; users who want a different behavior should edit the corresponding code in `main.jl`, or change the corresponding logging message to start with `[MAIN]` if a particular log message is wanted in multi-processing mode.

Users are encouraged to try out various combinations of numbers of **additional** processes (`-p`) and numbers of threads (`JULIA_NUM_THREADS`) for optimal throughput, and results are likely to differ for different computing systems.

Early tests on a 32-core node with 64G memory indicate that it ~16 processes with 2 threads each will likely still run, whereas more processes or threads will hit the memory limit soon for a regular configuration. Work to optimize the multi-process batch processing is ongoing. Note that applications built with RetrievalToolbox are generally quite memory-hungry, most objects persist to allow users to introspect most aspects of the retrieval at any given time.

For a given list of sounding IDs, either ingested via the `--sounding_id_list` argument, or by providing more than one ID through `--sounding_id`, the application chops up the list into equal lengths (if possible) and lets each process run through its own sub-list. Note that this is a static assignment for now, and due to the simple nature of this solution, there is no re-balancing of the workload **after** the application is called. So if one worker happens to be given a list of sounding IDs which all have a bad quality flag, that worker process will simply do nothing until all other processes are finished with their respective batch of retrievals. Also note that Julia's `SharedArray` functionality **only works on a single-node**, so running this set-up on a cluster where workers are spawned on different nodes will **not work**.

Should a worker process unexpectedly terminate, all references inside Julia's worker pool are removed automatically by Julia itself, so the program will keep running. However there is no re-balancing of the terminated worker's sounding ID list, so those soundings will not be processed in this run.

!!! note "This is proof-of-concept as of now!"
    The current implementation of parallel processing is mainly a proof-of-concept and has not been tested for large-scale production. Users who need a more sophisticated multi-processing set-ups need to implement their own solution.