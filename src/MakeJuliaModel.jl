"""

    make_julia_model(problem::VLJuliaModelObject; 
        intermediate_representation_dictionary::Union{Dict{String,Any},Nothing}=nothing,
        logger::Union{Nothing,SimpleLogger}=nothing) -> VLResult

This documentation is going to be awesome. Magical.
"""
function make_julia_model(problem::VLJuliaModelObject; 
    intermediate_representation_dictionary::Union{Dict{String,Any},Nothing}=nothing,
    logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

    # initialize -
    src_component_set = Set{NamedTuple}()
    network_component_set = Set{NamedTuple}()
    root_component_set = Set{NamedTuple}()
    
    try

        # check: do we have an intermediate representation?
        ir_dictionary = intermediate_representation_dictionary
        if (isnothing(intermediate_representation_dictionary) == true)
            
            # parse the vff document -
            result = read_model_document(problem)
            if (isa(result.value, Exception) == true)
                throw(result.value) # TODO: what is the diff between re-throw -vs- throw?
            end
            ir_dictionary = result.value
        end

        # Step 0: load the Defaults.toml file -
        path_to_model_file = problem.path_to_model_file
        path_to_defaults_file = joinpath(dirname(path_to_model_file), problem.defaults_file_name)
        default_dictionary = TOML.parsefile(path_to_defaults_file)
        ir_dictionary["default_dictionary"] = default_dictionary

        # Step 1: Generate the "distribution" files, copy to thier location -
        path_to_output_dir = problem.path_to_output_dir

        # build the problem dictionary -
        data_dictionary_component = generate_data_dictionary_program_component(problem, ir_dictionary)
        push!(src_component_set, data_dictionary_component)

        # build Kinetics.jl program component -
        kinetics_program_component = generate_kinetics_program_component(problem, ir_dictionary)
        push!(src_component_set, kinetics_program_component)

        # build Control.jl program component -
        control_program_component = generate_control_program_component(problem, ir_dictionary)
        push!(src_component_set, control_program_component)

        # build Balances.jl program component -
        balances_program_component = generate_balances_program_component(problem, ir_dictionary)
        push!(src_component_set, balances_program_component)

        # build the include program component -
        include_program_component = generate_include_program_component(problem, ir_dictionary)
        push!(root_component_set, include_program_component)

        # build the include program component -
        init_program_component = generate_init_program_component(problem, ir_dictionary)
        push!(root_component_set, init_program_component)

        # build the include program component -
        driver_program_component = generate_driver_program_component(problem, ir_dictionary)
        push!(root_component_set, driver_program_component)

        # build the stoichiometric_matrix -
        stoichiometric_matrix_component = generate_stochiometric_matrix_component(ir_dictionary)
        push!(network_component_set, stoichiometric_matrix_component)

        # build the degradation_matrix -
        degradation_matrix_component = generate_dilution_degradation_matrix_component(ir_dictionary)
        push!(network_component_set, degradation_matrix_component)

        # dump src components to disk -
        _output_path_to_src_distribution_files = joinpath(path_to_output_dir, "src")
        write_program_components_to_disk(_output_path_to_src_distribution_files, src_component_set)

        # dump root components to disk -
        write_program_components_to_disk(path_to_output_dir, root_component_set)

        # dump network components to disk -
        _output_path_to_network_distribution_files = joinpath(_output_path_to_src_distribution_files, "network")
        write_program_components_to_disk(_output_path_to_network_distribution_files, network_component_set)

        # return -
        return VLResult(nothing)
    catch error
        return VLResult(error)
    end
end