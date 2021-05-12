abstract type VLAbstractModelObject end

# concrete types -
struct VLResult{T}
    value::T
end

struct VLJuliaModelObject <: VLAbstractModelObject

    # data -
    path_to_model_file::String
    path_to_output_dir::String
    defaults_file_name::String 

    # constructor -
    function VLJuliaModelObject(path_to_model_file, path_to_output_dir; 
        defaults_file_path::String="Defaults.toml")
        
        # build new model object -
        this = new(path_to_model_file, path_to_output_dir, defaults_file_path)
    end
end

function check(result::VLResult; 
    logger::Union{Nothing,SimpleLogger}=nothing)::(Union{Nothing,T} where T <: Any)

    # ok, so check, do we have an error object?
    # Yes: log the error if we have a logger, then throw the error. 
    # No: return the result.value

    try

        # Error case -
        if (isa(result.value, Exception) == true)
            
            # get the error object -
            error_object = result.value

            # get the error message as a String -
            error_message = sprint(showerror, error_object)
        
            # log -
            if (isnothing(logger) == false)
                with_logger(logger) do
                    @error(error_message)
                end
            end

            # throw -
            throw(result.value)
        end

        # default -
        return result.value

    catch error
        rethrow(error)
    end
end