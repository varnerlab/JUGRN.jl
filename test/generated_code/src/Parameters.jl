# ----------------------------------------------------------------------------------- #
# Copyright (c) 2021 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the &quot;Software&quot;), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and&#x2F;or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #


function update_parameter_dictionary(parameter_dictionary::Dict{String,Any}, varargs...)::Union{Nothing, Dict{String,Any}}
    
    # Default: return the unmodified dictionary -
    @info "Default implementation of update_parameter_dictionary returns unmodified dictionary"

    # return parameter dictionary -
    return parameter_dictionary
end

function generate_default_parameter_dictionary()::Dict{String,Any}

    # initialize -
    parameter_dictionary = Dict{String,Any}()
    system_type_flag = 	:CF_TXTL

    try

        # load the stoichiometric_matrix (SM) and degradation_dilution_matrix (DM) -
        SM = readdlm(joinpath(_PATH_TO_NETWORK,"Network.dat"))
        DM = readdlm(joinpath(_PATH_TO_NETWORK,"Degradation.dat"))

        # build the species initial condition array -
        initial_condition_array = [
			5.0	;	#	1	gene_venus	units: nM
			0.0	;	#	2	mRNA_venus	units: nM
			0.0	;	#	3	P_venus	units: nM

			# translation capacity - 
			100.0	;	#	4	translation capacity	units: NA
		]

        # setup the system dimension -
        	
		number_of_transcription_processes = 1
		number_of_translation_processes = 1


        # build the system species concentration array -
        system_concentration_array = [
			0.07	;	#	RNAP	units: µM
			2.3	;	#	RIBOSOME	units: µM
			0.035	;	#	P_σ70	units: µM
		]

        # get the biophysical parameters for this system type -
        biophysical_parameters_dictionary = get_biophysical_parameters_dictionary_for_system(system_type_flag)

        # setup the model parameter array -
        model_parameter_array = [

			# control parameters: gene_venus
			0.001	;	#	1	W_gene_venus
			1.0		;	#	2	W_gene_venus_P_σ70
			1.0		;	#	3	K_gene_venus_P_σ70
			1.0		;	#	4	n_gene_venus_P_σ70

			# transcription kinetic limit parameters: gene_venus
			1.0		;	#	5	K_gene_venus
			1.0		;	#	6	𝛕_gene_venus

			# translation kinetic limit parameters: P_venus
			1.0		;	#	7	K_P_venus
			1.0		;	#	8	𝛕_P_venus

			# species degradation constants - 
			1.0		;	#	9	𝛳_mRNA_venus
			1.0		;	#	10	𝛳_P_venus
		]


        # setup the parameter symbol index map -
        model_parameter_symbol_index_map = Dict{Symbol,Int}()
		model_parameter_symbol_index_map[:W_gene_venus] = 1
		model_parameter_symbol_index_map[:W_gene_venus_P_σ70] = 2
		model_parameter_symbol_index_map[:K_gene_venus_P_σ70] = 3
		model_parameter_symbol_index_map[:n_gene_venus_P_σ70] = 4
		model_parameter_symbol_index_map[:K_gene_venus] = 5
		model_parameter_symbol_index_map[:𝛕_gene_venus] = 6
		model_parameter_symbol_index_map[:K_P_venus] = 7
		model_parameter_symbol_index_map[:𝛕_P_venus] = 8
		model_parameter_symbol_index_map[:𝛳_mRNA_venus] = 9
		model_parameter_symbol_index_map[:𝛳_P_venus] = 10


        # setup the inverse parameter symbol index map -
        inverse_model_parameter_symbol_index_map = Dict{Int, Symbol}()
		inverse_model_parameter_symbol_index_map[1] = :W_gene_venus
		inverse_model_parameter_symbol_index_map[2] = :W_gene_venus_P_σ70
		inverse_model_parameter_symbol_index_map[3] = :K_gene_venus_P_σ70
		inverse_model_parameter_symbol_index_map[4] = :n_gene_venus_P_σ70
		inverse_model_parameter_symbol_index_map[5] = :K_gene_venus
		inverse_model_parameter_symbol_index_map[6] = :𝛕_gene_venus
		inverse_model_parameter_symbol_index_map[7] = :K_P_venus
		inverse_model_parameter_symbol_index_map[8] = :𝛕_P_venus
		inverse_model_parameter_symbol_index_map[9] = :𝛳_mRNA_venus
		inverse_model_parameter_symbol_index_map[10] = :𝛳_P_venus


        # specific growth rate (default units: h^-1)
        μ = 0.0 # default units: h^-1

        # == DO NOT EDIT BELOW THIS LINE ========================================================== #
        parameter_dictionary["initial_condition_array"] = initial_condition_array
        parameter_dictionary["number_of_states"] = length(initial_condition_array)
        parameter_dictionary["number_of_transcription_processes"] = number_of_transcription_processes
        parameter_dictionary["number_of_translation_processes"] = number_of_translation_processes
        parameter_dictionary["system_concentration_array"] = system_concentration_array
        parameter_dictionary["biophysical_parameters_dictionary"] = biophysical_parameters_dictionary
        parameter_dictionary["model_parameter_array"] = model_parameter_array
        parameter_dictionary["model_parameter_symbol_index_map"] = model_parameter_symbol_index_map
        parameter_dictionary["inverse_model_parameter_symbol_index_map"] = inverse_model_parameter_symbol_index_map
        parameter_dictionary["stoichiometric_matrix"] = SM
        parameter_dictionary["dilution_degradation_matrix"] = DM
        parameter_dictionary["specific_growth_rate"] = μ

        # return -
        return parameter_dictionary
        # ========================================================================================= #
    catch error
        throw(error)
    end
end
