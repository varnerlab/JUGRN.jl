function _build_ic_array_snippet(model::VLJuliaModelObject, ir_dictionary::Dict{String,Any})::String

    # build the text fragments -
    ic_buffer=Array{String,1}()
    default_dictionary = ir_dictionary["default_dictionary"]

    # go ...
    +(ic_buffer,"initial_condition_array = [";suffix="\n")
    
    # iterate through species, build ic records -
    df_species = ir_dictionary["model_species_table"]
    (number_of_species,_) = size(df_species)
    for species_index = 1:number_of_species
        
        species_symbol = df_species[species_index,:symbol]
        species_type = df_species[species_index, :type]

        # what is the value?
        value = 0.0
        if (species_type == :DNA)
            value = default_dictionary["problem_dictionary"]["DNA_initial_condition"]
        elseif (species_type == :RNA)
            value = default_dictionary["problem_dictionary"]["RNA_initial_condition"]
        elseif (species_type == :PROTEIN)
            value = default_dictionary["problem_dictionary"]["PROTEIN_initial_condition"]
        end

        +(ic_buffer,"$(value)\t;\t#\t$(species_index)\t$(species_symbol)\tunits: nM\n"; prefix="\t\t\t")
    end

    +(ic_buffer,"]"; prefix="\t\t")
    
    # flatten and return -
    flat_ic_buffer = ""
    [flat_ic_buffer *= line for line in ic_buffer]
    return flat_ic_buffer
end

function _build_system_species_concentration_snippet(model::VLJuliaModelObject, ir_dictionary::Dict{String,Any})::String

    # initialize -
    system_buffer = Array{String,1}()

    # get the system block -
    system_dictionary_array = ir_dictionary["system_dictionary"]["list_of_system_species"]

    # here we go - open
    +(system_buffer, "system_concentration_array = ["; suffix="\n")
    for (index,system_species_dictionary) in enumerate(system_dictionary_array)
        
        # ok, block will be a dictionary with a value, and units keys - 
        value = system_species_dictionary["concentration"]
        units = system_species_dictionary["units"]
        symbol = system_species_dictionary["symbol"]

        # build the line -
        +(system_buffer, "$(value)\t;\t#\t$(symbol)\tunits: $(units)"; suffix="\n", prefix="\t\t\t")

    end

    # close -
    +(system_buffer,"]"; prefix="\t\t")


    # flatten and return -
    flat_buffer = ""
    [flat_buffer *= line for line in system_buffer]
    return flat_buffer
end

function _build_sequence_length_snippet(model::VLJuliaModelObject, ir_dictionary::Dict{String,Any})::String

    # initialize - 

end

function _build_model_species_alias_snippet(model::VLJuliaModelObject, ir_dictionary::Dict{String,Any})::String

    # initialize -
    buffer = Array{String,1}()

    # get the list of model species -
    df_species = ir_dictionary["model_species_table"]
    (number_of_species,_) = size(df_species)
    for species_index = 1:number_of_species
        
        # grab the symbol -
        species_symbol = df_species[species_index,:symbol]
        
        # this weird setup is because to the Mustache layout -
        if (species_index == 1)
            
            # build the record -
            +(buffer, "$(species_symbol) = x[$(species_index)]"; suffix="\n")
        elseif (species_index > 1 && species_index < number_of_species)
            
            # build the record -
            +(buffer, "$(species_symbol) = x[$(species_index)]"; prefix="\t", suffix="\n")
        else
            
            # build the record -
            +(buffer, "$(species_symbol) = x[$(species_index)]"; prefix="\t")
        end
    end

    # flatten and return -
    flat_buffer = ""
    [flat_buffer *= line for line in buffer]
    return flat_buffer
end

function _build_system_type_snippet(model::VLJuliaModelObject, ir_dictionary::Dict{String,Any})::String

    # initialize -
    buffer = Array{String,1}()

    # get the system block -
    system_type = ir_dictionary["system_dictionary"]["system_type"]

    # close -
    +(buffer,":$(system_type)"; prefix="\t")

    # flatten and return -
    flat_buffer = ""
    [flat_buffer *= line for line in buffer]
    return flat_buffer
end

# == MAIN METHODS BELOW HERE ======================================================================= #
function generate_data_dictionary_program_component(model::VLJuliaModelObject, ir_dictionary::Dict{String,Any})::NamedTuple

    # initialize -
    filename = "Problem.jl"
    buffer = Array{String,1}()
    template_dictionary = Dict{String,Any}()

    try 

        # build the snippets required by the 
        template_dictionary["copyright_header_text"] = build_julia_copyright_header_buffer(ir_dictionary)
        template_dictionary["initial_condition_array_block"] = _build_ic_array_snippet(model, ir_dictionary)
        template_dictionary["system_species_array_block"] = _build_system_species_concentration_snippet(model, ir_dictionary)
        template_dictionary["system_type_flag"] = _build_system_type_snippet(model, ir_dictionary)

        # write the template -
        template = mt"""
        {{copyright_header_text}}
        function generate_problem_dictionary()::Dict{String,Any}
            
            # initialize -
            problem_dictionary = Dict{String,Any}()
            system_type_flag = {{system_type_flag}}

            try

                # open a connection to the parameters db -


                # build the species initial condition array -
                {{initial_condition_array_block}}

                # build the system species concentration array -
                {{system_species_array_block}}

                


                # == DO NOT EDIT BELOW THIS LINE ======================================================= #
                problem_dictionary["initial_condition_array"] = initial_condition_array
                problem_dictionary["system_concentration_array"] = system_concentration_array

                # return -
                return problem_dictionary
                # ====================================================================================== #
            catch error
                throw(error)
            end
        end
        """

        # render step -
        flat_buffer = render(template, template_dictionary)
        
        # package up into a NamedTuple -
        program_component = (buffer=flat_buffer, filename=filename, component_type=:buffer)

        # return -
        return program_component
    catch error
        rethrow(error)
    end
end

function generate_kinetics_program_component(model::VLJuliaModelObject, ir_dictionary::Dict{String,Any})::NamedTuple

    # initialize -
    filename = "Kinetics.jl"
    template_dictionary = Dict{String,Any}()

    try

        # build snippets -
        template_dictionary["copyright_header_text"] = build_julia_copyright_header_buffer(ir_dictionary)
        template_dictionary["model_species_alias_block"] = _build_model_species_alias_snippet(model, ir_dictionary)

        # setup the template -
        template = mt"""
        {{copyright_header_text}}
        function calculate_transcription_kinetic_limit_array(t::Float64, x::Array{Float64,1}, 
            problem_dictionary::Dict{String,Any})::Array{Float64,1}
            
            # initialize -
            kinetic_limit_array = Array{Float64,1}()
            
            # alias the model species -
            {{model_species_alias_block}}

            
        
            # return -
            return kinetic_limit_array
        end

        function calculate_translation_kinetic_limit_array(t::Float64, x::Array{Float64,1}, 
            problem_dictionary::Dict{String,Any})::Array{Float64,1}
        
            # initialize -
            kinetic_limit_array = Array{Float64,1}()
            
            # alias the model species -
            {{model_species_alias_block}}

            
            
        
            # return -
            return kinetic_limit_array
        end

        function calculate_dilution_degradation_array(t::Float64, x::Array{Float64,1}, 
            problem_dictionary::Dict{String,Any})::Array{Float64,1}
            
            # initialize -
            degradation_dilution_array = Array{Float64,1}()
            
            # alias the model species -
            {{model_species_alias_block}}

            
            
        
            # return -
            return degradation_dilution_array
        end
        """

        # render step -
        flat_buffer = render(template, template_dictionary)
        
        # package up into a NamedTuple -
        program_component = (buffer=flat_buffer, filename=filename, component_type=:buffer)

        # return -
        return program_component
    catch error
        rethrow(error)
    end
end

function generate_balances_program_component(model::VLJuliaModelObject, 
    ir_dictionary::Dict{String,Any})::NamedTuple

    # initialize -
    filename = "Balances.jl"
    template_dictionary = Dict{String,Any}()

    try

        # build snippets -
        template_dictionary["copyright_header_text"] = build_julia_copyright_header_buffer(ir_dictionary)

        # setup the template -
        template = mt"""
        {{copyright_header_text}}
        function Balances(dx,x, problem_dictionary,t)

            # get 
            number_of_states = problem_dictionary["number_of_states"]
            AM = problem_dictionary["dilution_degradation_matrix"]
            SM = problem_dictionary["stoichiometric_matrix"]
             
            # calculate the TX and TL kinetic limit array -
            transcription_kinetic_limit_array = calculate_transcription_kinetic_limit_array(t,x,problem_dictionary)
            translation_kinetic_limit_array = calculate_translation_kinetic_limit_array(t,x,problem_dictionary)
            
            # calculate the TX and TL control array -
            u = calculate_transcription_control_array(t,x,problem_dictionary)
            w = calculate_translation_control_array(t,x,problem_dictionary)

            # calculate the rate of transcription and translation -
            r_TX = transcription_kinetic_limit_array.*u
            r_TL = translation_kinetic_limit_array.*w
            rV = [r_TX ; r_TL]

            # calculate the degradation and dilution rates -
            r_dd = calculate_dilution_degradation_array(t,x,problem_dictionary)

            # compute the model equations -
            dxdt = SM*rV + AM*r_dd

            # package -
            for index = 1:number_of_states
                dx[index] = dxdt[index]
            end
        end
        """

        # render step -
        flat_buffer = render(template, template_dictionary)
        
        # package up into a NamedTuple -
        program_component = (buffer=flat_buffer, filename=filename, component_type=:buffer)

        # return -
        return program_component
    catch error
        rethrow(error)
    end
end
# ================================================================================================= #