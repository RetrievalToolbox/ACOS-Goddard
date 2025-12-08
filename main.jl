# The retrieval toolkit - every function you use with a RE.xxx is from there
# Any other function is either from a third-party module or from this repo.

using RetrievalToolbox
const RE = RetrievalToolbox

using ArgParse
using Dates
using Distributed # needed for batch processing with shared memory
using DocStringExtensions
using HDF5
using Interpolations
using LinearAlgebra
using Logging, LoggingExtras
using LoopVectorization
using Printf
using ProgressMeter
using Statistics
using StatsBase
using Unitful
using TimerOutputs


# Functions to read in data
include("io.jl")
# Function to create a state vector
include("state_vector.jl")
# Function that processes a scenen through the inversion algorithm
include("process_scene.jl")
# Forward model
include("forward_model.jl")
# Non-uniform sampling
include("NUS.jl")
# Routines to process the results
include("collect_results.jl")

logger = ConsoleLogger(stderr, Logging.Info);
global_logger(logger);

# During multi-processing, let us reduce the number of log messages printed out
function multi_filter(log_args)
    startswith(log_args.message, "[MAIN]")
end

if nprocs() > 1
    global_logger(ActiveFilteredLogger(multi_filter, global_logger()));
end


"""
Splits a vector into near-equal chunks. If the vector is shorter than the number of
workers, then empty sub-lists will be returned.
"""
function distribute_work(arr, nsub = nprocs())

    chunk_size = div(length(arr), nsub)
    remainder = length(arr) % nsub

    return [arr[i*chunk_size + 1 + min(i, remainder):(i+1)*chunk_size + min(i+1, remainder)]
            for i in 0:nsub-1]
end

#=
    Global variables
=#

# Number of pressure levels on the retrieval grid
N_RT_lev = 20
# What number type we use (keep this Float64)
#=
    Note - for flexibility, RetrievalToolbox does not use many eplicit Float64 variables
    or arrays. Users could thus, in theory, write a retrieval application that makes use
    of mostly Float32 to save memory. XRTM, however, explicitly requires Float64s in many
    functions, hence we have no choice but to leave this as Float64 for the time being.
=#
my_type = Float64

# This is where the cmdline arguments are processed
include("args.jl")


function main(barrier_channel, sync_channel, ARGS_in)

    #=
        Prepare the paths that contain the input data
    =#

    # Parse the command line arguments into a dictionary
    args = parse_commandline(ARGS_in)

    # Check if we have a sounding ID list AND a single sounding ID
    if !isempty(args["sounding_id"]) & !isnothing(args["sounding_id_list"])
        @error "Cannot have both `sounding_id` and `sounding_id_list`!"
        exit(1)
    elseif isempty(args["sounding_id"]) & isnothing(args["sounding_id_list"])
        @error "Must have either `sounding_id` or `sounding_id_list`!"
        exit(1)
    end

    if !isempty(args["sounding_id"])
        sounding_id_list = args["sounding_id"] # these are already integers
    else
        # Parse the sounding ID list
        sounding_id_list_txt = readlines(args["sounding_id_list"])
        # This needs to be parsed into integers
        sounding_id_list = parse.(Ref(Int), sounding_id_list_txt)
    end

    # Each worker grabs its own share:
    sounding_id_list = distribute_work(sounding_id_list)[myid()]

    @info "[MAIN] ACOS-Goddard will process N=$(length(sounding_id_list)) scenes."

    # Parse the spectral windows
    spec_array = [parse(Int, x) for x in split(args["spec"], ",")]

    # Which spectrometers are we processing?
    if 1 in spec_array
        @info "Processing O2 A-band"
    end
    if 2 in spec_array
        @info "Processing Weak CO2 band"
    end
    if 3 in spec_array
        @info "Processing Strong CO2 band"
    end

    #=
        #########################
        Create the needed objects
        #########################
        The following objects are very "read-only" from the viewpoint of the rest
        of the algorithm, hence they can be loaded once into memory *before* the sounding
        ID loop.
    =#

    # Read in spectroscopy, depending on the window configuration. We
    # keep them in a Dictionary that are accessed via a simple string
    # so that e.g. `abscos["CO2"]` gets you the CO2 absco.

    #=
        This is the ONLY part in the code, for now, in which we are making use of Julia's
        distributed computing function. The ABSCO tables are a large chunk of the total
        memory usage in this application. We thus let only the root process read them into
        memory, and then share it to all workers. This is possible, because the
        ABSCOSpectroscopy4D objects all accept SharedArrays for the coefficient array.
        This way, we do not duplicate the array in memory, but ALL workers use the same
        memory space.
    =#

    local absco_channel

    if myid() == 1
        # Root creates the absco coordination channel
        absco_channel = RemoteChannel(() ->
            Channel{Dict{String, RE.AbstractSpectroscopy}}(1))
        # Signal to all workers that channel is ready
        put!(barrier_channel, absco_channel)
    else
        # Worker takes it from the channel
        absco_channel = take!(barrier_channel)
        # Puts it back for the next one to grab it
        put!(barrier_channel, absco_channel)
    end


    if myid() == 1

        # Generate SharedArray only if we have multi-processing..
        distributed = nprocs() > 1

        # Root channel creates the ABSCO dict, and loads the data..
        abscos = Dict{String, RE.AbstractSpectroscopy}()

        # Find out which ABSCO loader to use.. (this is from io.jl)
        load_spectroscopy = which_ABSCO_loader(args["o2_spec"])

        # If we use the O2 A-band, we need O2 at least
        if 1 in spec_array
            @info "Reading in O2 spectroscopy ..."
            abscos["O2"] = load_spectroscopy(
                # Pass the path to the ABSCO file
                args["o2_spec"],
                spectral_unit=:Wavelength,
                distributed=distributed
            )
            @info "... done!"

            # Apply a user-defined spectroscopy scaling factor. This can be done to the
            # entire spectroscopic table since oxygen is only present in the A-band
            # anyway.
            # Note: scaling should only be done by ONE worker in case of distributed
            @info "Scaling O2 cross sections by $(args["o2_scale"])"
            abscos["O2"].cross_section[:] .*= args["o2_scale"]

        end

        # If we use either of the two CO2 bands, we need H2O and CO2
        if (2 in spec_array) | (3 in spec_array)
            @info "Reading in CO2 spectroscopy ..."

            # Find out which ABSCO loader to use.. (this is from io.jl)
            load_spectroscopy = which_ABSCO_loader(args["co2_spec"])

            abscos["CO2"] = load_spectroscopy(
                # Pass the path to the ABSCO file
                args["co2_spec"],
                spectral_unit=:Wavelength,
                distributed=distributed
            )
            @info "... done!"

            # Apply a user-defined scale factor for CO2 spectroscopy in the weak band
            @info "Scaling CO2 cross sections for weak CO2 by " *
                "$(args["co2_scale_weak"])"
            idx_weak = findall(abscos["CO2"].ww * abscos["CO2"].ww_unit .< 2.0u"µm")
            abscos["CO2"].cross_section[idx_weak,:,:,:] .*= args["co2_scale_weak"]


            # Apply a user-defined scale factor for CO2 spectroscopy in the strong band
            @info "Scaling CO2 cross sections for strong CO2 by " *
                "$(args["co2_scale_strong"])"
            idx_strong = findall(abscos["CO2"].ww * abscos["CO2"].ww_unit .> 2.0u"µm")
            abscos["CO2"].cross_section[idx_strong,:,:,:] .*= args["co2_scale_strong"]


            @info "Reading in H2O spectroscopy ..."

            # Find out which ABSCO loader to use.. (this is from io.jl)
            load_spectroscopy = which_ABSCO_loader(args["h2o_spec"])

            abscos["H2O"] = load_spectroscopy(
                # Pass the path to the ABSCO file
                args["h2o_spec"],
                spectral_unit=:Wavelength,
                distributed=distributed
            )
            @info "... done!"
        end

        # In case we are running with more than one process
        if !isempty(workers())
            # Put the ABSCO dict into the remote channel
            @info "[MAIN] Root puts ABSCO into channel.."
            put!(absco_channel, abscos)
        end
    else
        # Other workers wait their turn and recieve them
        abscos = take!(absco_channel)
        @info "[MAIN] .. received ABSCO dictonary!"
        # And put it back for the next one to receive..
        put!(absco_channel, abscos)
        @info "[MAIN] .. puts ABSCO dictonary back into channel."

    end # End myid() == 1

    # This section is only needed for more than one process:
    if !isempty(workers())
        # Let all workers catch up
        if myid() > 1
            @info "[MAIN] .. sends sync signal"
            put!(sync_channel, myid())
        else
            slist = Int[]
            for id in 2:nprocs()
                @info "[MAIN] Waits for signal from $(id)"
                push!(slist, take!(sync_channel))
                @info "[MAIN] Recieved sync signal from $(id)"
            end
        end
    end

    # Distributed:
    # From here on, ALL processes have a valid reference to ABSCO, but thanks to the
    # SharedArray, the coefficient arrays only exist once in memory.

    #=
        For every ABSCO object, we create a new gas object. Same as the ABSCOs,
        we stick them into a dictionary that is accessed by the gas name, e.g.
        `gases["O2"]` delivers the oxygen gas object which itself uses the
        `abscos["O2"]` ABSCO table.
    =#

    gas_units = Dict(
        "O2" => Unitful.NoUnits,
        "H2O" => Unitful.percent,
        "CO2" => Unitful.ppm
        )

    gases = Dict{String, RE.GasAbsorber}()

    for (gas_name, gas_absco) in abscos
        @info "Creating gas $(gas_name) with ABSCO $(gas_absco.file_name) ..."

        gases[gas_name] = RE.GasAbsorber(
            gas_name, # Name of the gas
            gas_absco, # ABSCO object
            zeros(my_type, N_RT_lev), # Placeholder for VMR level profile
            gas_units[gas_name]
        )

        @info "... done!"
    end


    @info "Reading in solar model from $(args["solar_model"]) .. "

    if args["solar_model"][end-2:end] == ".h5"
        # Use the ACOS solar model
        solar_models = Dict(
                spec => RE.OCOHDFSolarModel(
                # Path to the solar model file
                args["solar_model"],
                # "spec" referring to OCO-HDF solar model band index
                spec,
                spectral_unit=:Wavelength
            ) for spec in spec_array
        )

    elseif args["solar_model"][end-2:end] == ".nc"
        # Use the TSIS model
        solar_models = Dict(
            spec => RE.TSISSolarModel(
                args["solar_model"],
                spectral_unit=:Wavelength
            ) for spec in spec_array
        )

        for (spec, sm) in solar_models
            # In-place conversion of solar model units
            RE.convert_solar_model_to_photons!(sm)
        end
    end
    @info "... done!"


    #=
        Create the spectral windows!
        (they are created using a helper function will will make sure
        that the spectral grid is congruent with the ABSCO sampling)

        These window limits could be read from a config file as well..
    =#

    spectral_windows = Dict{Int, RE.SpectralWindow}()

    if 1 in spec_array

        spectral_windows[1] = RE.spectralwindow_from_ABSCO(
            "o2",
            0.7592152960425568, 0.7714610963322822, 0.760,
            5.0e-3,
            abscos["O2"],
            Unitful.μm # or u"μm"
        )

    end

    if 2 in spec_array
        spectral_windows[2] = RE.spectralwindow_from_ABSCO(
            "weak_co2",
            1.5979699810546206, 1.617757907862207, 1.60,
            5.0e-3,
            abscos["CO2"],
            Unitful.μm # or u"μm"
        )
    end

    if 3 in spec_array
        spectral_windows[3] = RE.spectralwindow_from_ABSCO(
            "strong_co2",
            2.0476103795182246, 2.0796702096182256, 2.0,
            5.0e-3,
            abscos["CO2"],
            Unitful.μm # or u"μm"
        )
    end

    # Need some static data (NUS grid, CO2 covariance matrix)
    h5_oco_static = h5open("./example_data/l2_oco_static_input.h5", "r")

    #=
        Non-uniform sampling (NUS) grid
        (read the NUS knots from the OCO static file)
    =#

    nus_dict = Dict{Int, Vector{Bool}}()
    for spec in spec_array
        @info "Creating non-uniform sampling data for $(spectral_windows[spec])."
        # The NUS knots are in wavenumber
        nus_ww = h5_oco_static["Spectrum_Sampling/nonuniform_grid_$(spec)"][:]
        # .. convert to wavelength
        nus_wl = 1e4 ./ nus_ww
        # find which wavelengths in this spectral window are in common with the
        # NUS knots from the file
        wl_common = intersect(nus_wl, spectral_windows[spec].ww_grid)
        # Where are the common wavelengths found?
        nus_idx = searchsortedfirst.(Ref(spectral_windows[spec].ww_grid), wl_common);

        # Now create a new boolean mask array with the same length as the
        # high-res wavelengths
        nus = zeros(Bool, spectral_windows[spec].N_hires);
        # .. and set it to `true` for every spectral point that is to be calculated!
        nus[nus_idx] .= 1;
        # Store it in the dict.
        nus_dict[spec] = nus;
    end


    # Grab the CO2 covariance matrix from the ACOS static file
    # (and multiply by (ppm^2) to get to ppm units.
    CO2_covar = h5_oco_static["Gas/CO2/covariance"][:,:] .* 1e12

    close(h5_oco_static)

    # Big sounding ID loop! Note this is wrapped in a @sync as we need all workers to
    # finish at the same time. Otherwise the sharedarray from the root process would
    # cease to exist while others might still be running..

    for sounding_id in sounding_id_list

        # Free up memory, maybe useful for multi-processing
        # NOTE This is generally not recommended for fast retrievals, but ACOS retrievals
        # take ~minutes.
        GC.gc()

        @info "##################################"
        @info "[MAIN] RETRIEVING $(sounding_id)"
        @info "##################################"

        # Touch an init file if needed
        if args["touch_init"]
            touch_fname = joinpath(args["output"], "$(sounding_id).h5.init")
            run(`touch $(touch_fname)`)
        end
        #=
            Read in the retrieval inputs for this scene
        =#

        @info "Reading in scene inputs ..."

        scene_inputs = Dict(
            spec => OCO_read_inputs_from_l1b_and_met(
                args["L1b"],
                args["L2Met"],
                band_idx=spec,
                sounding_id=sounding_id,
                l2cpr_fname=args["L2CPr"]
                ) for spec in spec_array
            )

        @info "... done!"

        # Skip bad soundings immediately
        _input = first(values(scene_inputs))
        if _input["sounding_qual_flag"][sounding_id] != 0
            @info "[MAIN] Bad sounding quality for $(sounding_id). Exiting."

            # Remove init file
            if args["touch_init"]
                touch_fname = joinpath(args["output"], "$(sounding_id).h5.init")
                new_touch_fname = joinpath(args["output"], "$(sounding_id).h5.error")
                run(`mv $(touch_fname) $(new_touch_fname)`)
            end

            # Move on to next sounding
            continue
        end

        #=
            Create dispersion objects!
        =#
        my_dispersions = Dict{Int, RE.SimplePolynomialDispersion}()
        for spec in spec_array
            @info "Creating dispersion object for $(spectral_windows[spec])."

            # Note that the dispersion polynomial coefficients are present
            # in every instance of `scene_inputs`.

            my_dispersions[spec] = RE.SimplePolynomialDispersion(
                scene_inputs[spec]["dispersion"][spec],
                1:1016,
                spectral_windows[spec]
            );
        end


        #=
            Create the atmospheric elements to be used in the retrieval!
        =#

        # Start with an empty vector
        atm_elements = RE.AbstractAtmosphereElement[]

        # Add all available gases
        for (gas_name, gas) in gases
            push!(atm_elements, gas)
        end

        # Add Rayleigh scattering
        push!(atm_elements, RE.RayleighScattering())

        #=
            Add the aerosols we use for every scene (ice cloud, water cloud, statospheric
            sulphur). This can change from scene to scene, so we must do it for each
            new sounding ID.
        =#

        if args["aerosols"]

            @info "Using aerosols."

            #=
                Ice cloud
            =#
            @info "Adding ice cloud ..."
            miemom_ice = read_ACOS_aerosol(
                "example_data/l2_aerosol_combined.h5",
                "ice_cloud_MODIS6_deltaM_50",
                wavelength=true
            )
            aer_ice = RE.GaussAerosol(
                "Ice",
                miemom_ice,
                true, # relative pressure?
                0.0, # height p/psurf
                0.04, # width p/psurf
                Unitful.NoUnits,
                0.755, # ref wavelength
                u"µm", # ref wavelength unit
                exp(-5.382), # Total AOD at reference point
                );

            push!(atm_elements, aer_ice)

            #=
                Water cloud
            =#
            @info "Adding water cloud ..."
            miemom_water = read_ACOS_aerosol(
                "example_data/l2_aerosol_combined.h5",
                "wc_008",
                wavelength=true
            )
            aer_water = RE.GaussAerosol(
                "Water",
                miemom_water,
                true, # relative pressure?
                0.75, # height p/psurf
                0.1, # width p/psurf
                Unitful.NoUnits,
                0.755, # ref wavelength
                u"µm", # ref wavelength unit
                exp(-4.382), # Total AOD at reference point
                );

            push!(atm_elements, aer_water)

            #=
                Stratospheric aerosol
            =#
            @info "Adding stratospheric aerosol ..."
            miemom_strat = read_ACOS_aerosol(
                "example_data/l2_aerosol_combined.h5",
                "strat",
                wavelength=true
                );
            aer_strat = RE.GaussAerosol(
                "ST",
                miemom_strat,
                true, # relative pressure?
                0.03, # height p/psurf
                0.04, # width p/psurf
                Unitful.NoUnits,
                0.755, # ref wavelength
                u"µm", # ref wavelength unit
                exp(-5.116), # Total AOD at reference point
            );

            push!(atm_elements, aer_strat)

            #=
                Now add the two largest contributors to tropospheric aerosols from MERRA
            =#
            # These are potentially hard-coded, but really look at "/Metadata/CompositeAerosolTypes"
            trop_types = ["DU", "SS", "BC", "OC", "SO"]

            # (add 1 to account for 0-based vs 1-based indexing..)
            trop_idx = scene_inputs[first(spec_array)]["met"][sounding_id]["aerosol_sort_met"] .+ 1
            # Find where the indices are 1 or 2 (to get the top two aerosols)
            # (or adjust to your liking, ..)
            trop_choice_idx = indexin([1,2], trop_idx)
            trop_chosen = trop_types[trop_choice_idx]

            @info "Top two aerosols are: $(trop_chosen)"

            for idx in trop_choice_idx
                @info "Reading aerosol properties for $(trop_types[idx])"

                this_trop_miemom = read_ACOS_aerosol(
                    "example_data/l2_aerosol_combined.h5",
                    trop_types[idx],
                    wavelength=true
                )

                # Extract height (fraction of psurf), width (fraction of psurf), AOD
                gauss_params = scene_inputs[first(spec_array)]["met"][sounding_id]["aerosol_gauss_params_met"][:, idx]

                # The position of the Gaussian parameters was inferred from the L2FP code..
                # (might need checking). Note we make conversions here to make sure that all
                # inputs to the GaussAerosol have the correct float type.

                this_trop_height = convert(my_type, gauss_params[2])
                this_trop_width = convert(my_type, gauss_params[3])
                # Override -- ACOS says we have to do 0.05 width
                this_trop_width = convert(my_type, 0.05)
                this_trop_aod = convert(my_type, gauss_params[4])

                @info "Creating Gauss-type aerosol for $(trop_types[idx])."
                this_trop_gauss = RE.GaussAerosol(
                    trop_types[idx],
                    this_trop_miemom,
                    true, # relative pressure?
                    this_trop_height, # height p/psurf
                    this_trop_width, # width p/psurf
                    Unitful.NoUnits,
                    convert(my_type, 0.755), # ref wavelength
                    u"µm", # ref wavelength unit
                    this_trop_aod, # Total AOD at reference point
                )
                # Push this aerosol into the list of atmospheric elements

                push!(atm_elements, this_trop_gauss)

            end
        end
        # Finished populating the list of atmospheric elements!
        #@info "Atmospheric elements: $(atm_elements)"

        #@info "Creating state vector ..."
        # Create a state vector with a helper function.
        state_vector = create_statevector(
            spectral_windows,
            my_dispersions,
            atm_elements;
            lambertian=args["Lambertian"],
            retrieve_aerosols=args["retrieve_aerosols"],
            retrieve_psurf=args["retrieve_psurf"]
        )

        N_sv = length(state_vector)

        #=

            Create buffers. Buffers are pre-allocated objects and arrays which the various
            functions then use to in-place modify and place intermediate results into.
            This dramatically reduces the amount of allocations which would negatively
            impact overall performance.

        =#

        @info "Create buffers ..."

        # Will contain outputs of ISRF application. These values seem to work for
        # a normal 3-band OCO-type retrieval.
        inst_buf = RE.InstrumentBuffer(
            zeros(my_type, 5000),
            zeros(my_type, 5000),
            zeros(my_type, 1016),
        )
        # Will contain outputs of radiance, jacobians and
        # dispersion indices. We make this a ScalarRTBuffer because
        # the instrument only measures Intensity, and not the polarization components.

        # We use ScalarRTBuffer because the instrument measures scalar radiance!
        rt_buf = RE.ScalarRTBuffer(
            Dict(spectral_windows[i] => my_dispersions[i] for i in spec_array),
            RE.ScalarRadiance(my_type, 3*1016), # Hold the radiance
            Dict(sve => RE.ScalarRadiance(my_type, 3*1016) for sve in state_vector.state_vector_elements),
            Dict(spectral_windows[i] => zeros(Int, 0) for i in spec_array), # Hold the detector indices
            u"ph/s/m^2/sr/µm", # Radiance units for the forward model, this must be in the same units as the L1B data..
        )


        # Choose radiance type depending on whether user wants polarization or not
        if args["polarized"]
            rad_type = RE.VectorRadiance
        else
            rad_type = RE.ScalarRadiance
        end


        # Choose the surface type depending on user input
        if args["Lambertian"]
            # This is NOT the default, the number represents the polynomial order maximum
            surface_list = [(:Lambertian, 3) for x in spec_array]
        else
            #=
            This is the default as per ACOS specification. The first numer is the polynomial
            order maximum. The other 3 parameters are the RPV kernel parameters:
                1) Hot-spot parameter ρ
                2) Asymmetry parameter Θ
                3) Anisotropy parameter k

            Note that this order is kept from the way how XRTM ingests them.
            =#
            surface_list = [(:RPV, 3, 0.05, -0.1, 0.75) for x in spec_array]
        end


        # Will contain everything needed to run a retrieval
        buf = RE.EarthAtmosphereBuffer(
            state_vector, # The state vector
            [spectral_windows[i] for i in spec_array], # The spectral window (or a list of them)
            surface_list, # Surface type for each band (same order as spectral windows)
            atm_elements, # Atmospheric elements
            Dict(spectral_windows[i] => solar_models[i] for i in spec_array), # The solar model(s),
            [:XRTM for i in spec_array], # RT model
            rad_type, # Radiance type, depending on user choice
            rt_buf, # The pre-allocated RT buffer
            inst_buf, # The pre-allocated instrument buffer
            N_RT_lev, # The number of retrieval or RT pressure levels
            72, # The number of meteorological pressure levels
            my_type # The chosen Float data type (e.g. Float16, Float32, Float64)
        )

        @info "... done!"


        # Setting the RT options

        # Grab the number of high streams from ARGS
        Nhigh = args["Nhigh"]

        # The monochromatic calculations are doing using first-order polarized RT, and then
        # a fast two-stream multiple scattering scalar calculation.
        mo1 = Dict()
        mo2 = Dict()

        mo1["solvers"] = ["single"]
        mo1["add"] = true
        mo1["streams"] = 2
        mo1["options"] = [
            "output_at_levels",
            "calc_derivs",
            "source_solar",
            "psa",
            "sfi",
            "save_leg_polys",
            "save_phase_mats"
            ]

        # Make SS contributions polarized if needed
        if args["polarized"]
            push!(mo1["options"], "vector")
        end

        if args["aerosols"]
            push!(mo1["options"], "delta_m")
            push!(mo1["options"], "n_t_tms")
        end

        mo2["solvers"] = ["two_stream"]
        mo2["add"] = true
        mo2["streams"] = 2
        mo2["options"] = [
            "output_at_levels",
            "calc_derivs",
            "source_solar",
            "psa",
            "sfi",
            "save_leg_polys",
            "save_phase_mats"
            ]

        if args["aerosols"]
            push!(mo2["options"], "delta_m")
            push!(mo2["options"], "n_t_tms")
        end

        for (swin, rt) in buf.rt

            empty!(rt.model_options)
            push!(rt.model_options, mo1)
            push!(rt.model_options, mo2)

        end


        #=
            High-stream options for LSI
        =#

        hmo1 = Dict()
        hmo1["solvers"] = ["single"]
        hmo1["add"] = true
        hmo1["streams"] = Nhigh
        hmo1["options"] = [
            "output_at_levels",
            "calc_derivs",
            "source_solar",
            "psa",
            "sfi",
            "save_leg_polys",
            "save_phase_mats"
            ]

        # Make SS contributions polarized if needed
        if args["polarized"]
            push!(hmo1["options"], "vector")
        end

        if args["aerosols"]
            push!(hmo1["options"], "delta_m")
            push!(hmo1["options"], "n_t_tms")
        end

        hmo2 = Dict()
        hmo2["solvers"] = ["two_os"]
        hmo2["add"] = true
        hmo2["streams"] = Nhigh
        hmo2["options"] = [
            "output_at_levels",
            "calc_derivs",
            "source_solar",
            "vector",
            "psa",
            "sfi",
        ]

        hmo3 = Dict()
        hmo3["solvers"] = ["eig_bvp"]
        hmo3["add"] = true
        hmo3["streams"] = Nhigh
        hmo3["options"] = [
            "output_at_levels",
            "calc_derivs",
            "source_solar",
            "psa",
            "sfi",
            ]

        if args["aerosols"]
            push!(hmo3["options"], "delta_m")
            push!(hmo3["options"], "n_t_tms")
        end

        if args["LSI"]

            if args["polarized"]
                # Use SS (vector), 2OS (vector) and EIG_BVP (scalar) for polarized run
                high_options = [hmo1, hmo2, hmo3]
            else
                # Use SS (scalar) and EIG_BVP (scalar) for scalar run
                high_options = [hmo1, hmo3]
            end

        else

            # No LSI, set this to `nothing`
            high_options=nothing

        end

        #=
            Process the scene
        =#

        (solver, fm_kwargs), proc_time = @timed process_snid(
            spectral_windows,
            scene_inputs,
            state_vector,
            sounding_id,
            buf;
            max_iter=args["max_iterations"],
            gamma=args["gamma"],
            dsigma_scale=args["dsigma_scale"],
            co2_cov=CO2_covar,
            high_options=high_options,
            nus_dict=nus_dict,
            )

        @info "[MAIN] $(sounding_id) finished in $(proc_time) seconds"

        # Return buffer, solver and forward model arguments for single-ID retrieval
        if (args["max_iterations"] == 0) & (length(sounding_id_list) == 1)
            @info "No iterations requested - returning buffer, solver and FM keyword args."
            return buf, solver, fm_kwargs
        end

        this_result = collect_results(solver, buf)

        if !isnothing(this_result)

            # Output file name:
            outfname = joinpath(args["output"], "$(sounding_id).h5")

            results = Dict()
            results[sounding_id] = this_result
            @info "Writing results for $(sounding_id)"
            @debug "Opening $(outfname) for write access."
            fid = h5open(outfname , "w")
            for (k,v) in this_result
                isnothing(v) && continue

                @debug "Writing out $(k)"
                try
                    write(fid, k, v)
                catch
                    @warn "Failed writing $(k)."
                end
            end
            close(fid)

            # Remove init file
            if args["touch_init"]
                touch_fname = joinpath(args["output"], "$(sounding_id).h5.init")
                run(`rm $(touch_fname)`)
            end

            @info "[MAIN] Done retrieving $(sounding_id)"

        end

        # Return buffer, solver and forward model arguments for single-ID retrieval
        if (length(sounding_id_list) == 1)
            @info "Returning!"
            return buf, solver, fm_kwargs
        end

    end # End sounding ID loop

end

