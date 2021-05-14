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

        +(ic_buffer,"$(value)\t;\t#\t$(species_index)\t$(species_symbol)\n"; prefix="\t\t\t")
    end

    +(ic_buffer,"]"; prefix="\t\t",suffix="\n")
    
    flat_ic_buffer = ""
    [flat_ic_buffer *= line for line in ic_buffer]
    return flat_ic_buffer
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

        # write the template -
        template = mt"""
        {{copyright_header_text}}
        function generate_problem_dictionary()::Dict{String,Any}
            
            # initialize -
            problem_dictionary = Dict{String,Any}()

            try

                # build the initial condition array -
                {{initial_condition_array_block}}

                # == DO NOT EDIT BELOW THIS LINE ======================================================= #
                problem_dictionary["initial_condition_array"] = initial_condition_array
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
        throw(error)
    end
end
# ================================================================================================= #