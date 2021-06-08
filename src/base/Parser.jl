function _build_species_table(model_dictionary::Dict{String,Any})::DataFrame

    # initialize -
    species_table = DataFrame(symbol=String[], type=Symbol[], compartment=String[], sequence=Union{Missing,String,FASTA.Record}[])

    try

        # ok, so lets build an intermediate representation that is a bunch of DataFrames -
        list_of_species_dictionaries = model_dictionary["list_of_model_species"]
        for species_dictionary in list_of_species_dictionaries
            
            # grab -
            local_data_row = Union{String,Symbol,Missing,FASTA.Record}[]
            push!(local_data_row, species_dictionary["symbol"])
            push!(local_data_row, Symbol(species_dictionary["type"]))
            push!(local_data_row, species_dictionary["compartment"])

            # Lets load the sequence -
            seq_path = species_dictionary["sequence"]
            if (seq_path == "")
                push!(local_data_row, missing)
            else 
                open(FASTA.Reader, seq_path) do reader
                    for record in reader
                        push!(local_data_row, record)
                    end
                end
            end
            
            # push into data frame -
            push!(species_table, tuple(local_data_row...))
        end

        # return -
        return species_table
    catch error
        rethrow(error)
    end
end

# == PUBLIC METHODS BELOW HERE ==================================================================== #

"""
    read_model_document(model::VLJuliaModelObject) -> VLResult

This documentation is going to be awesome. 
"""
function read_model_document(model::VLJuliaModelObject)::VLResult

    # initialize -
    intermediate_representation_dictionary = Dict{String,Any}()

    try

        # where is the main file located?
        path_to_model_file = model.path_to_model_file

        # load the model JSON file -
        json_model_dictionary = JSON.parsefile(path_to_model_file)

        # get model species, transscription and translation models -
        model_species_table = _build_species_table(json_model_dictionary)
        list_of_transcription_models = json_model_dictionary["list_of_transcription_models"]
        list_of_translation_models = json_model_dictionary["list_of_translation_models"]
        system_type = json_model_dictionary["system"]["system_type"]
        
        # how many species do we have?
        (number_of_species, _) = size(model_species_table)

        # setup up the system dimension -
        intermediate_representation_dictionary["number_of_species"] = number_of_species
        intermediate_representation_dictionary["number_of_transcription_models"] = length(list_of_transcription_models)
        intermediate_representation_dictionary["number_of_translation_models"] = length(list_of_translation_models)
        intermediate_representation_dictionary["system_type"] = system_type

        # Build a species table -
        intermediate_representation_dictionary["model_species_table"] = model_species_table

        # Grab the system dictionary -
        intermediate_representation_dictionary["system_dictionary"] = json_model_dictionary["system"]

        # Grab the transcription models -
        intermediate_representation_dictionary["list_of_transcription_models"] = list_of_transcription_models

        # Grab the translation models -
        intermediate_representation_dictionary["list_of_translation_models"] = list_of_translation_models

        # return -
        return VLResult(intermediate_representation_dictionary)
    catch error
        return VLResult{Exception}(error)
    end
end
# == PUBLIC METHODS ABOVE HERE ==================================================================== #