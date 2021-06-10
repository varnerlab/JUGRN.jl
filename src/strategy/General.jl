function generate_stochiometric_matrix_component(ir_dictionary::Dict{String,Any})::NamedTuple

    # initialize -
    number_of_species = ir_dictionary["number_of_species"]
    system_type = ir_dictionary["system_type"]
    number_of_transcription_models = ir_dictionary["number_of_transcription_models"]
    number_of_translation_models = ir_dictionary["number_of_translation_models"]
    number_of_rates = number_of_transcription_models + number_of_translation_models
    
    # get the species table -
    model_species_table = ir_dictionary["model_species_table"]

    try
        
        # what is the overall dimension of the system -
        (number_of_species, _) = size(model_species_table)

        # how many genes do we have?
        number_of_genes = filter(x -> x == :DNA, model_species_table[:,:type]) |> length
        
        # blocks -
        gene_block = zeros(number_of_genes, number_of_rates)
        id_block = Matrix{Float64}(I, (number_of_species - number_of_genes), (number_of_species - number_of_genes))
        
        # stm -
        stoichiometric_matrix = [gene_block ; id_block]

        # ok, if we have a CF_ system type, then we have an extra row of zeros -
        if (contains(system_type,"CF_") == true)
            zero_row = zeros(1,number_of_rates)
            stoichiometric_matrix = [stoichiometric_matrix ; zero_row]
        end
        
        # package -
        stm_program_component = (matrix = stoichiometric_matrix, filename = "Network.dat", component_type = :matrix)

        # return -
        return stm_program_component
    catch error
        rethrow(error)
    end
end

function generate_dilution_degradation_matrix_component(ir_dictionary::Dict{String,Any})::NamedTuple
    
    # initialize -
    number_of_species = ir_dictionary["number_of_species"]
    system_type = ir_dictionary["system_type"]
    number_of_transcription_models = ir_dictionary["number_of_transcription_models"]
    number_of_translation_models = ir_dictionary["number_of_translation_models"]
    number_of_rates = number_of_transcription_models + number_of_translation_models
    
    # get the species table -
    model_species_table = ir_dictionary["model_species_table"]

    try 

        # what is the overall dimension of the system -
        (number_of_species, _) = size(model_species_table)

        # how many genes do we have?
        number_of_genes = filter(x -> x == :DNA, model_species_table[:,:type]) |> length
        
        # blocks -
        gene_block = zeros(number_of_genes, number_of_rates)
        id_block = -1 * Matrix{Float64}(I, (number_of_species - number_of_genes), (number_of_species - number_of_genes))
        
        # stm -
        degradation_matrix = [gene_block ; id_block]

        # ok, if we have a CF_ system type, then we have an extra row of zeros -
        if (contains(system_type,"CF_") == true)
            zero_row = zeros(1,number_of_rates)
            degradation_matrix = [degradation_matrix ; zero_row]
        end
        
        # package -
        degradation_program_component = (matrix = degradation_matrix, filename = "Degradation.dat", component_type = :matrix)

        # return -
        return degradation_program_component
    catch error
        rethrow(error)
    end
end