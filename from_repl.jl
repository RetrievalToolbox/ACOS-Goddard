#=
This code snippet demonstrates how one can supply command line arguments into the
application while running from a Julia REPL. There must be a list of arguments,
`my_args` here, in which the arguments and values are listed as strings. The order of
argument-value pairs does not matter, however it is required that the corresponding
value must follow the argument descriptor, just like in a normal command line.

The list of strings is then added to the command line argument variable of the currently
running Julia session.

=#

using Distributed

my_args = [
    #"--solar_model", "./data/hybrid_reference_spectrum_p005nm_resolution_c2022-11-30_with_unc.nc",
    "--solar_model", "./example_data/l2_solar_model.h5", # Path to solar model file
    "--L1b", "./example_data/2021030111564431_inputs.h5", # Path to L1B file
    "--L2Met", "./example_data/2021030111564431_inputs.h5", # Path to L2Met file
    "--L2CPr", "./example_data/2021030111564431_inputs.h5", # Path to L2CPr file
    "--sounding_id", "2021030111564431", # Which single sounding ID to process
    #"--sounding_id", "2021030111564431", "2021030111564431", # Which single sounding ID(s) to process
    #"--sounding_id_list", "sounding_id_list.txt", # Which soundings (list) to process
    "--output", "./", # Path to where the output file will be written
    "--touch_init", "true", # Whether to write init-files
    "--o2_spec", "./data/o2_v52.hdf", # path to O2 ABSCO file
    "--o2_scale", "1.0048", # scaling factor for oxygen spectroscopy
    "--co2_spec", "./data/co2_v52.hdf", # path to CO2 ABSCO file
    "--co2_scale_weak", "0.994", # scaling factor for CO2 spectroscopy in weak band
    "--co2_scale_strong", "0.998", # scaling factor for CO2 spectroscopy in strong band
    "--h2o_spec", "./data/h2o_v52.hdf", # path to H2O ABSCO file
    "--spec", "1,2,3", # spectrometers to use
    "--polarized", "true", # use polarization in radiative transfer?
    "--Lambertian", "false", # use Lambertian surface instead of RPV BRDF?
    "--aerosols", "true", # add aerosols to the atmosphere?
    "--retrieve_aerosols", "true", # retrieve aerosols?
    "--retrieve_psurf", "true", # retrieve surface pressure (only works if A-band is retrieved)?
    "--LSI", "true", # use LSI to get high-accuracy MS calculations?
    "--Nhigh", "16", # what number of (half-)streams to use with LSI?
    "--gamma", "500.0", # LM-gamma initial value
    "--dsigma_scale", "2.0", # scale factor for convergence check (higher means faster convergence)
    "--max_iterations", "10", # number of iterations before the retrieval is forced to stop
]

# For distributed, we need to have channels to synchronize operations and moving the
# spectroscopy data around:
barrier_channel = RemoteChannel(() -> Channel{RemoteChannel}(1))
sync_channel = RemoteChannel(() -> Channel{Any}(nworkers() - 1))

if nprocs() > 1
    # For multi-processing we must @everywhere all these calls
    @everywhere my_args = $my_args
    @everywhere empty!(ARGS)
    @everywhere push!(ARGS, my_args...)

    @everywhere include("main.jl")
    @everywhere main($barrier_channel, $sync_channel, ARGS)
else
    # In case of single processing, users might want to inspect the results contents of
    # the buffer, the solver and the forward model keyword arguments.
    empty!(ARGS)
    push!(ARGS, my_args...)

    include("main.jl")
    main(barrier_channel, sync_channel, ARGS)
end