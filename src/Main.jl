function make(problem::VLJuliaModelObject; 
    intermediate_representation_dictionary::Union{Dict{String,Any},Nothing}=nothing,
    logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

    # initialize -
    # ...
    
    try

        # check: do we have an intermediate representation?
        ir_dictionary = intermediate_representation_dictionary
        if (isnothing(intermediate_representation_dictionary) == true)
            
            # parse the vff document -
            result = parse_model_document(problem)
            if (isa(result.value, Exception) == true)
                rethrow(result.value) # TODO: what is the diff between re-throw -vs- throw?
            end
            ir_dictionary = result.value
        end

        # Step 1: Copy the "distribution" files to thier location -
        path_to_output_dir = problem.path_to_output_dir
        
        


        # return -
        return VLResult(nothing)
    catch error
        return VLResult(error)
    end
end