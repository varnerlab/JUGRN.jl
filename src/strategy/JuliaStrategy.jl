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

function _build_system_species_alias_snippet(model::VLJuliaModelObject, ir_dictionary::Dict{String,Any})::String

    # initialize -
    buffer = Array{String,1}()

    # get the system_dictionary -
    system_dictionary = ir_dictionary["system_dictionary"]    
    list_of_system_species = system_dictionary["list_of_system_species"]
    for (index,species_dictionary) in enumerate(list_of_system_species)
        
        # get data -
        symbol = species_dictionary["symbol"]

        # formulate system species line -
        line = "$(symbol) = system_array[$(index)]"

        # add to buffer -
        +(buffer, line; prefix="\t", suffix="\n")
    end

    



    # flatten and return -
    flat_buffer = ""
    [flat_buffer *= line for line in buffer]
    return flat_buffer
end

function _build_u_variable_snippet(model::VLJuliaModelObject, ir_dictionary::Dict{String,Any})::String

    # initialize -
    buffer = Array{String,1}()

    # get the list of transcription models -
    list_of_transcription_models = ir_dictionary["list_of_transcription_models"]
    number_of_transcription_models = length(list_of_transcription_models)
    for (index,model_dictionary) in enumerate(list_of_transcription_models)
        
        # get data from the model dictionary -
        title = model_dictionary["title"]
        input = model_dictionary["input"]
        output = model_dictionary["output"]
        polymerase_symbol = model_dictionary["polymerase_symbol"]

        # comment string -
        comment_string = "# $(title): $(input) -[$(polymerase_symbol)]- $(output)"

        # for this promoyer - get the list of activators
        initialize_activator_set = "$(title)_activator_set = Array{Float64,1}()\n"
        initialize_activator_set *= "\tpush!($(title)_activator_set, 1.0)"

        list_of_activators = model_dictionary["list_of_activators"]
        local_activator_buffer = ""
        for (index, activator_dictionary) in enumerate(list_of_activators)
            
            # get activator symbol -
            activator_symbol = activator_dictionary["symbol"]
            activator_type = activator_dictionary["type"]

            # check: if sigma factor, then actor is bound to RNAP
            if (activator_type == "positive_sigma_factor")
                local_activator_buffer *= "$(activator_symbol)_$(polymerase_symbol) = $(polymerase_symbol)*f($(activator_symbol),K,n)\n"
                local_activator_buffer *= "\tpush!($(title)_activator_set, $(activator_symbol)_$(polymerase_symbol))"
            else
                # nothibng for now ...
            end
        end
        numerator_string = "A = sum(W.*$(title)_activator_set)"

        # for this promoter - process the list of repressors 
        list_of_repressors = model_dictionary["list_of_repressors"]
        initialize_repressor_set = "$(title)_repressor_set = Array{Float64,1}()"
        local_repressor_buffer = ""
        for (index, repressor_dictionary) in enumerate(list_of_repressors)
        
            # get repressor symbol -
            repressor_symbol = repressor_dictionary["symbol"]
            repressor_type = repressor_dictionary["type"]

            # check the type -
            if (repressor_type == "negative_protein_repressor")
                
                # if we get here, then we have a counter agent -
                counter_agent = repressor_dictionary["counter_agent"]
                
                # build the activate component -
                local_repressor_buffer *= "$(repressor_symbol)_active = $(repressor_symbol)*(1.0 - f($(counter_agent), K, n))\n"
                local_repressor_buffer *= "\tpush!($(title)_repressor_set, $(repressor_symbol)_active)"

            else
                # nothing for now -
            end
        end

        repressor_term_list = "R = 0"
        if (isempty(list_of_repressors) == false)
            repressor_term_list = "R = sum(W.*$(title)_repressor_set)"
        end 


        # layout -
        if (index==1)
        
            # setup comment -
            +(buffer, comment_string; suffix="\n")
        
        else
            
            # setup comment -
            +(buffer, comment_string; suffix="\n", prefix="\t")
        end

        
        # write the activator lines -
        if (isempty(list_of_activators) == false)
            +(buffer, "# $(title) activator set"; prefix="\t", suffix="\n")
            +(buffer, initialize_activator_set; prefix="\t", suffix="\n")
            +(buffer, local_activator_buffer; prefix="\t", suffix="\n")
            +(buffer, numerator_string; prefix="\t", suffix="\n")
            
            if (isempty(list_of_repressors) == false)
                +(buffer,"\n")
            end
        end


        if (isempty(list_of_repressors) == false)
            +(buffer, "# $(title) repressor set"; prefix="\t", suffix="\n")
            +(buffer, initialize_repressor_set; prefix="\t", suffix="\n")
            +(buffer, local_repressor_buffer; prefix="\t", suffix="\n")
            +(buffer, repressor_term_list; prefix="\t", suffix="\n")  
        else 
            +(buffer, repressor_term_list; prefix="\t", suffix="\n") 
        end

        +(buffer, "push!(u_array, u(A,R))"; suffix="\n", prefix="\t")
    end

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

            # system dimensions and structural matricies -
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
            rX = transcription_kinetic_limit_array.*u
            rL = translation_kinetic_limit_array.*w
            rV = [rX ; rL]

            # calculate the degradation and dilution rates -
            rd = calculate_dilution_degradation_array(t,x,problem_dictionary)

            # compute the model equations -
            dxdt = SM*rV + AM*rd

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

function generate_control_program_component(model::VLJuliaModelObject, 
    ir_dictionary::Dict{String,Any})::NamedTuple

    # initialize -
    filename = "Control.jl"
    template_dictionary = Dict{String,Any}()

    try

        # build snippets -
        template_dictionary["copyright_header_text"] = build_julia_copyright_header_buffer(ir_dictionary)
        template_dictionary["u_variable_snippet"] = _build_u_variable_snippet(model,ir_dictionary)
        template_dictionary["model_species_alias_block"] = _build_model_species_alias_snippet(model, ir_dictionary)
        template_dictionary["system_species_alias_block"] = _build_system_species_alias_snippet(model, ir_dictionary)

        # build template -
        template=mt"""
        {{copyright_header_text}}
        
        # binding function -
        f(x,K,n) = (x^n)/(K^n+x^n)
        u(A,R) = A/(1+A+R)
        
        # calculate the u-variables -
        function calculate_transcription_control_array(t::Float64, x::Array{Float64,1}, 
            problem_dictionary::Dict{String,Any})::Array{Float64,1}

            # initialize -
            number_of_transcription_processes = problem_dictionary["number_of_transcription_processes"]
            u_array = Array{Float64,1}(undef,number_of_transcription_processes)

            # alias the state vector -
            {{model_species_alias_block}}

            # alias the system species -
            system_array = problem_dictionary["system_concentration_array"]
            {{system_species_alias_block}}

            # == CONTROL LOGIC BELOW ================================================================= #
            {{u_variable_snippet}}
            # == CONTROL LOGIC ABOVE ================================================================= #

            # return -
            return u_array
        end

        # calculate the w-variables -
        function calculate_translation_control_array(t::Float64, x::Array{Float64,1}, 
            problem_dictionary::Dict{String,Any})::Array{Float64,1}
            
            # defualt: w = 1 for all translation processes -
            number_of_translation_processes = problem_dictionary["number_of_translation_processes"]
            w = ones(number_of_translation_processes)
            
            # return -
            return w
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

function generate_include_program_component(model::VLJuliaModelObject, ir_dictionary::Dict{String,Any})::NamedTuple

    # initialize -
    filename = "Include.jl"
    template_dictionary = Dict{String,Any}()

    try

        # build snippets -
        template_dictionary["copyright_header_text"] = build_julia_copyright_header_buffer(ir_dictionary)

        # build template -
        template=mt"""
        {{copyright_header_text}}

        # setup paths -
        const _PATH_TO_ROOT = pwd()
        const _PATH_TO_SRC = joinpath(_PATH_TO_ROOT, "src")

        # get packages -
        import Pkg
        Pkg.activate(_PATH_TO_ROOT)
        Pkg.add(name="DifferentialEquations")
        Pkg.add(name="VLModelParametersDB")
        
        # use packages -
        using DifferentialEquations
        using VLModelParametersDB

        # include my codes -
        include("./src/Problem.jl")
        include("./src/Balances.jl")
        include("./src/Kinetics.jl")
        include("./src/Control.jl")
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