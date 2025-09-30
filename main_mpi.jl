#=
    NOTE that this script is no longer supported and should be phased out!
=#


# The retrieval toolkit - every function you use with a RE.xxx is from there
# Any other function is either from a third-party module or from this repo.
using RetrievalToolbox
const RE = RetrievalToolbox


using ArgParse
using CSV
using DataFrames
using Dates
using DocStringExtensions
using HDF5
using Interpolations
using LinearAlgebra
using Logging, LoggingExtras
using LoopVectorization
using MPI
using Printf
using ProgressMeter
using Statistics
using StatsBase
using UnicodePlots
using Unitful

#=
Initialize MPI and determine which Rank this process is.
=#

MPI.Init()
MPI_comm = MPI.COMM_WORLD

MPI_rank = MPI.Comm_rank(MPI_comm)
MPI_size = MPI.Comm_size(MPI_comm)

# MPI functions
include("mpi_helpers.jl")
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



#=
    Global variables
=#

# Number of pressure levels on the retrieval grid
N_RT_lev = 20
# what number type we use (keep this Float64)
my_type = Float64

# This is where the cmdline arguments are processed
include("args.jl")


function main()

    #=
        Prepare the paths that contain the input data
    =#

    # Parse the command line arguments into a dictionary
    args = parse_commandline()

    #=
    We define a new logger here, which conveniently keeps track of the Rank of this process.
    =#

    my_logger(logger) = TransformerLogger(logger) do log
        merge(log, (; message = "($(MPI_rank)) $(log.message)"))
    end

    if isnothing(args["logfile"])
        Logging.ConsoleLogger(stdout, Logging.Info) |> my_logger |> global_logger
    else
        logger_io = open(args["logfile"], "w+")
        Logging.ConsoleLogger(logger_io, Logging.Info) |> my_logger |> global_logger
    end

    if isnothing(args["sounding_id_list"])
        @error "Must have the sounding ID list argument!"
    end

    # Parse the spectral windows
    spec_array = [parse(Int, x) for x in split(args["spec"], ",")]

    # Determine which gases we have, based on the spectral window configuration
    used_gases = []
    if 1 in spec_array
        push!(used_gases, "O2")
    end
    if (2 in spec_array) | (3 in spec_array)
        push!(used_gases, "CO2")
        push!(used_gases, "H2O")
    end

    #=
        Retrieval inputs are held in dictionaries, which we initiliaze here, so that they
        are present for all ranks.
    =#
    abscos = Dict{String, Any}() # gas name -> absco object
    gases = Dict{String, Any}() # gas name -> gas object
    spectral_windows = Dict{Int, Any}() # band index -> spectral window object
    nus_dict = Dict{Int, Any}() # band index -> NUS dictionary
    solar_models = Dict{Int, Any}() # band index -> solar model data
    oco_static = Dict{String, Any}() # data path -> anything
    aerosol_data = Dict{String, Any}() # data path -> anything

    #=
        Set-up and initial IO will be done by rank = 0 only
    =#

    if MPI_rank == 0

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
            Read in the retrieval inputs
        =#

        @info "Reading in scene inputs ..."

        # For now it seems to be simply quicker to read in the entire L1B file, rather
        # than picking specific sounding IDs.
        all_scene_inputs = Dict(
            spec => OCO_read_inputs_from_l1b_and_met(
                args["L1b"],
                args["L2Met"],
                band_idx=spec,
                l2cpr_fname=args["L2CPr"]
                ) for spec in spec_array
            )

        @info "... done!"


        #=
                Here we check whether we are given a list of sounding IDs, via the
                argument "--sounding_id_list", or a single sounding ID, via the argument
                "--sounding_id".
        =#
        @info "Reading list of sounding IDs .."
        sounding_id_list = parse.(Ref(Int), readlines(args["sounding_id_list"]))
        N_ids = length(sounding_id_list)
        @info "About to process N = $(N_ids) soundings."

        # Split the sounding ID list into a roughly equal number of chunks to be processed
        # by each active rank.
        if MPI_size > 1
            #slicer_bounds = collect(0:Int(round(N_ids / (MPI_size-1))):N_ids)
            slicer_bounds = convert.(Ref(Int),
                round.(collect(range(0, N_ids, length=MPI_size+1)))
            )
        else
            slicer_bounds = [0, N_ids]
        end

        if slicer_bounds[end] != N_ids
            slicer_bounds[end] = N_ids
        end

        for i in 0:MPI_size - 1
            @info "Rank $(i): $(slicer_bounds[i+1]+1) through $(slicer_bounds[i+2])"
        end

        # How many soundings are in the L1B file?
        N_ids_infile = length(keys(all_scene_inputs[spec_array[1]]["fp_frame_idx"]))

        # Split up the scene inputs into dictionaries that can be sent to ranks
        send_scene_inputs = Dict{Any, Any}()
        # .. these will be send_scene_inputs[rank][band][...]
        for rank in 0:MPI_size - 1

            send_scene_inputs[rank] = Dict{Any, Any}()
            s_start = slicer_bounds[rank+1] + 1
            s_end = slicer_bounds[rank+2]

            for spec in spec_array

                send_scene_inputs[rank][spec] = Dict{Any, Any}()

                # Loop through all parent keys
                for pkey in keys(all_scene_inputs[spec])
                    send_scene_inputs[rank][spec][pkey] = Dict{Any, Any}()

                    # Let's identify whether this is a sounding ID or not
                    if length(keys(all_scene_inputs[spec][pkey])) == N_ids_infile
                        for snid in sounding_id_list[s_start:s_end]
                            send_scene_inputs[rank][spec][pkey][snid] =
                                all_scene_inputs[spec][pkey][snid]
                        end
                        #@info "Key $(pkey) has sounding IDs!"
                    else
                        # This is probably a footprint-related quantity, so grab all
                        send_scene_inputs[rank][spec][pkey] = all_scene_inputs[spec][pkey]
                        #@info "Key $(pkey) has no sounding information!"
                    end
                end

            end
        end

        # No need to retain "all_scene_inputs" anymore, let the GC free this space..
        all_scene_inputs = nothing

        #=
            #########################
            Create the needed objects
            #########################
        =#

        # For the OCO-type solar model
        @info "Reading in solar model from $(args["solar_model"]) .. "

        for spec in spec_array

            if args["solar_model"][end-2:end] == ".h5"
                solar_models[spec] = RE.OCOHDFSolarModel(
                        # Path to the solar model file
                        args["solar_model"],
                        # "spec" referring to OCO-HDF solar model band index
                        spec,
                        spectral_unit=:Wavelength
                    )

            elseif args["solar_model"][end-2:end] == ".nc"

                solar_models[spec] = RE.TSISSolarModel(
                    args["solar_model"],
                    spectral_unit=:Wavelength
                )

                RE.convert_solar_model_to_photons!(solar_models[spec])

            end

        end

        @info "... done!"


        # Read in spectroscopy, depending on the window configuration.

        # If we use the O2 A-band, we need O2 at least
        if 1 in spec_array
            @info "Reading in O2 spectroscopy ..."
            abscos["O2"] = RE.load_ABSCO_spectroscopy(
                # Pass the path to the ABSCO file
                "./data/o2_v5.1_v180916+mlawer181031.hdf",
                spectral_unit=:Wavelength
            )
            @info "... done!"

            # Apply a user-defined spectroscopy scaling factor. This can be done to the
            # entire spectroscopic table since oxygen is only present in the A-band
            # anyway.
            @info "Scaling O2 cross sections by $(args["o2_scale"])"
            @turbo abscos["O2"].cross_section[:] .*= args["o2_scale"]

        end

        # If we use either of the two CO2 bands, we need H2O and CO2
        if (2 in spec_array) | (3 in spec_array)
            @info "Reading in CO2 spectroscopy ..."
            abscos["CO2"] = RE.load_ABSCO_spectroscopy(
                # Pass the path to the ABSCO file
                "./data/co2_v5.2_test.hdf",
                spectral_unit=:Wavelength
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
            abscos["H2O"] = RE.load_ABSCO_spectroscopy(
                # Pass the path to the ABSCO file
                "./data/h2o_atm18+mtckd32.hdf",
                spectral_unit=:Wavelength
            )
            @info "... done!"
        end


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

        # Read full file into a dict
        h5_oco_static = h5open("./example_data/l2_oco_static_input.h5", "r")
        walk_h5!(h5_oco_static, oco_static)
        close(h5_oco_static)

        h5_aerosol = h5open("./example_data/l2_aerosol_combined.h5", "r")
        walk_h5!(h5_aerosol, aerosol_data)
        close(h5_aerosol)

        #=
            Non-uniform sampling (NUS) grid
            (read the NUS knots from the OCO static file)
        =#

    end

    # Let MPI catch up
    MPI.Barrier(MPI_comm)

    #=
        Below, we broadcast all common objects needed for the retrieval.
    =#



    #=
        Broadcast gases
    =#
    if MPI_rank > 0

        # Create empty placeholder objects
        for gas_name in used_gases
            gases[gas_name] = nothing
        end
    end

    for (gas_name, gas) in gases
        MPI_rank == 0 && @info "Broadcasting gas $(gas_name)"
        gases[gas_name] = MPI.bcast(gas, MPI_comm; root=0)
        #gases[gas_name] = MPI.Bcast(gas, 0, MPI_comm)
    end

    #=
        Broadcast solar model
    =#
    if MPI_rank > 0

        # Create empty placeholder objects
        for spec in spec_array
            solar_models[spec] = nothing
        end
    end

    MPI_rank == 0 && @info "Broadcasting solar models"
    for spec in spec_array
        solar_models[spec] = MPI.bcast(solar_models[spec], MPI_comm; root=0)
    end

    #=
        Broadcast static data, and aerosols
    =#
    MPI_rank == 0 && @info "Broadcasting ACOS static data"
    oco_static = MPI.bcast(oco_static, MPI_comm; root=0)
    MPI_rank == 0 && @info "Broadcasting ACOS aerosol data"
    aerosol_data = MPI.bcast(aerosol_data, MPI_comm; root=0)

    # Let MPI catch up
    MPI.Barrier(MPI_comm)


    # Send scene data to ranks 1 .. MPI_size-1
    if MPI_rank == 0

        for dest in 1:MPI_size-1
            tag = MPI_tags[:scene_inputs][dest]
            @info "Sending scene inputs to $(dest) via tag $(tag)"
            sreq = MPI.isend(send_scene_inputs[dest], MPI_comm; dest=dest, tag=tag)
            wait(sreq)
        end

        # Retain for rank 0
        scene_inputs = send_scene_inputs[0]


    # receive data from rank 0
    else
        tag = MPI_tags[:scene_inputs][MPI_rank]
        scene_inputs = MPI.recv(MPI_comm; source=0, tag=tag)
        @info "Recieved scene inputs via tag $(tag)"
    end

    MPI.Barrier(MPI_comm)


    #=
        Create the spectral windows!
        (they are created using a helper function will will make sure
        that the spectral grid is congruent with the ABSCO sampling)
    =#
    if 1 in spec_array

        spectral_windows[1] = RE.spectralwindow_from_ABSCO(
            "o2",
            0.7592152960425568, 0.7714610963322822, 0.760,
            5.0e-3,
            gases["O2"].spectroscopy,
            Unitful.μm # or u"μm"
        )

    end

    if 2 in spec_array
        spectral_windows[2] = RE.spectralwindow_from_ABSCO(
            "weak_co2",
            1.5979699810546206, 1.617757907862207, 1.60,
            5.0e-3,
            gases["CO2"].spectroscopy,
            Unitful.μm # or u"μm"
        )
    end

    if 3 in spec_array
        spectral_windows[3] = RE.spectralwindow_from_ABSCO(
            "strong_co2",
            2.0476103795182246, 2.0796702096182256, 2.0,
            5.0e-3,
            gases["CO2"].spectroscopy,
            Unitful.μm # or u"μm"
        )
    end


    for spec in spec_array
        @info "Creating non-uniform sampling data for $(spectral_windows[spec])."
        # The NUS knots are in wavenumber
        nus_ww = oco_static["/Spectrum_Sampling/nonuniform_grid_$(spec)"][:]
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
    # (and divide by (ppm^2) to get to ppm units.
    CO2_covar = oco_static["/Gas/CO2/covariance"][:,:] .* 1e12


    #=
        Create dispersion objects!
    =#
    my_dispersions = Dict{Int, Any}()
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
        ##############
        Sounding loop!
        ##############

        For ACOS retrievals, each scene can have a different state vector, due to the fact
        that the ingested aerosol types are scene-dependent. Since the retrievals take
        minutes, rather than fractions of a second, we are taking a performance hit here
        by creating the buffer object for each scene. If the algorithm ever becomes fast
        enough for this to be a significant bottleneck, other options should be explored.

    =#

    # Which sounding IDs are we processing? Grab those from a sounding-id dependent
    # dictionary..
    spec_first = first(spec_array)
    snid_process_list = collect(keys(scene_inputs[spec_first]["met"]))
    sort!(snid_process_list)

    for sounding_id in snid_process_list

        # This will be the output name
        out_fname = "$(args["output"])/$(sounding_id).h5"

        if isfile("$(out_fname).init") | isfile("$(out_fname).error") | isfile(out_fname)
            @info "Skipping sounding $(sounding_id) (exists)"
            continue
        else
            run(`touch $(out_fname).init`)
        end

        @info "Processing sounding $(sounding_id)"


        if scene_inputs[spec_first]["sounding_qual_flag"][sounding_id] != 0
            @info "Bad sounding quality. Skipping."
            run(`mv $(out_fname).init $(out_fname).bad`)
            continue
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
            sulphur)
        =#

        if args["aerosols"]

            @info "Using aerosols."
            #=
                Ice cloud
            =#
            @info "Adding ice cloud ..."
            miemom_ice = read_ACOS_aerosol(
                "data/l2_aerosol_combined.h5",
                "ice_cloud_MODIS6_deltaM_50",
                wavelength=true
            )
            aer_ice = RE.GaussAerosol(
                "Ice",
                miemom_ice,
                true,
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
                "data/l2_aerosol_combined.h5",
                "wc_008",
                wavelength=true
            )
            aer_water = RE.GaussAerosol(
                "Water",
                miemom_water,
                true,
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
                "data/l2_aerosol_combined.h5",
                "strat",
                wavelength=true
                );
            aer_strat = RE.GaussAerosol(
                "ST",
                miemom_strat,
                true,
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
            # These are potentially hard-coded, but really look at
            # "/Metadata/CompositeAerosolTypes"
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
                    "data/l2_aerosol_combined.h5",
                    trop_types[idx],
                    wavelength=true
                )

                # Extract height (fraction of psurf), width (fraction of psurf), AOD
                gauss_params = scene_inputs[first(spec_array)]["met"][sounding_id]["aerosol_gauss_params_met"][:, idx]

                # The position of the Gaussian parameters was inferred from the L2FP code..
                # (might need checking)
                this_trop_height = convert(my_type, gauss_params[2])
                this_trop_width = convert(my_type, gauss_params[3])
                # Override -- ACOS says we have to do 0.05 width
                this_trop_width = 0.05
                this_trop_aod = convert(my_type, gauss_params[4])

                @info "Creating Gauss-type aerosol for $(trop_types[idx])."
                this_trop_gauss = RE.GaussAerosol(
                    trop_types[idx],
                    this_trop_miemom,
                    true,
                    this_trop_height, # height p/psurf
                    this_trop_width, # width p/psurf
                    Unitful.NoUnits,
                    0.755, # ref wavelength
                    u"µm", # ref wavelength unit
                    this_trop_aod, # Total AOD at reference point
                )
                # Push this aerosol into the list of atmospheric elements

                push!(atm_elements, this_trop_gauss)

            end
        end
        # Finished populating the list of atmospheric elements!
        @info "Atmospheric elements: $(atm_elements)"


        @info "Creating state vector ..."
        # Create a state vector with a helper function.
        state_vector = create_statevector(
            spectral_windows,
            my_dispersions,
            atm_elements;
            retrieve_aerosols=args["retrieve_aerosols"]
        )

        N_sv = length(state_vector)

        #=

            Create buffers.
            Buffers are pre-allocated objects and arrays which the various
            functions then use to in-place modify and place intermediate
            results into. This dramatically reduces the amount of allocations
            which would negatively impact overall performance.

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
            rad_type, # radiance type, depending on user choice
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

            We wrap this portion into a try/catch so that the MPI process
            does not crash if something goes wrong.
        =#
        try

            global solver, fm_kwargs = process_snid(
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

            if isnothing(solver)
                @warn "Failed performing retrieval for $(sounding_id)"
                run(`rm $(out_fname).init`)
                run(`touch $(out_fname).error`)
                continue
            end

        catch err
            # Signal error and move on
            @warn "Failed performing retrieval for $(sounding_id)"
            @warn exception=(err, catch_backtrace())
            run(`rm $(out_fname).init`)
            run(`touch $(out_fname).error`)
            continue
        end

        if args["max_iterations"] > 0

            results = Dict()
            try
                results[sounding_id] = collect_results(solver, buf)
            catch
                @warn "Collecting results failed!"
                run(`rm $(out_fname).init`)
                run(`touch $(out_fname).error`)
                continue
            end
            failed = false

            @info "Opening $(out_fname) for write access."
            fid = h5open(out_fname, "w")
            for (k,v) in results[sounding_id]
                isnothing(v) && continue
                @info "Writing out $(k)"
                try
                    write(fid, k, v)
                catch
                    @warn "Failed writing $(k)."
                    # Signal error and move on
                    failed = true
                    break
                end
            end
            close(fid)

            if failed
                @warn "Failed writing outputs."
                run(`rm $(out_fname).init`)
                run(`touch $(out_fname).error`)
                continue
            end


            @info "Done."

            # Get rid of init file
            run(`rm $(out_fname).init`)
        else
            @time solver.forward_model(solver.state_vector)
        end

    end # End sounding ID loop!
    ###########################

    # Let MPI catch up
    MPI.Barrier(MPI_comm)


end
