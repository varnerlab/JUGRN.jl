function read_model_document(model::VLJuliaModelObject)::VLResult

    # initialize -
    intermediate_representation_dictionary = Dict{String,Any}()

    try

        

        # return -
        return VLResult(intermediate_representation_dictionary)
    catch error
        return VLResult{Exception}(error)
    end
end