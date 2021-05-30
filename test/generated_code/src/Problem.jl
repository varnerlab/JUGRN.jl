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

function generate_problem_dictionary()::Dict{String,Any}

    # initialize -
    problem_dictionary = Dict{String,Any}()
    system_type_flag = 	:CF_TXTL

    try

        # load the stoichiometric_matrix (SM) and degradation_dilution_matrix (DM) -
        SM = readdlm(joinpath(_PATH_TO_NETWORK,"Network.dat"))
        DM = readdlm(joinpath(_PATH_TO_NETWORK,"Degradation.dat"))

        # build the species initial condition array -
        initial_condition_array = [
			5.0	;	#	1	gene_gntR	units: nM
			5.0	;	#	2	gene_venus	units: nM
			0.0	;	#	3	mRNA_gntR	units: nM
			0.0	;	#	4	mRNA_venus	units: nM
			0.0	;	#	5	P_gntR	units: nM
			0.0	;	#	6	P_venus	units: nM
		]

        # setup the system dimension -
        	
		number_of_transcription_processes = 2
		number_of_translation_processes = 2


        # build the system species concentration array -
        system_concentration_array = [
			0.07	;	#	RNAP	units: µM
			0.07	;	#	RIBOSOME	units: µM
			1.0	;	#	P_σ70	units: µM
			1.0	;	#	M_gluconate_c	units: µM
		]

        # get the biophysical parameters for this system type -
        biophysical_parameters_dictionary = get_biophysical_parameters_dictionary_for_system(system_type_flag)

        # setup the model parameter array -
        model_parameter_array = [

			# control parameters: gene_venus
			0.001	;	#	1	W_gene_venus
			1.0		;	#	2	W_gene_venus_P_σ70
			1.0		;	#	3	W_gene_venus_P_gntR
			1.0		;	#	4	K_gene_venus_P_σ70
			1.0		;	#	5	n_gene_venus_P_σ70
			1.0		;	#	6	K_gene_venus_P_gntR
			1.0		;	#	7	n_gene_venus_P_gntR

			# control parameters: gene_gntR
			0.001	;	#	8	W_gene_gntR
			1.0		;	#	9	W_gene_gntR_P_σ70
			1.0		;	#	10	K_gene_gntR_P_σ70
			1.0		;	#	11	n_gene_gntR_P_σ70

			# transcription kinetic limit parameters: gene_venus
			1.0		;	#	12	K_gene_venus
			1.0		;	#	13	𝛕_gene_venus

			# transcription kinetic limit parameters: gene_gntR
			1.0		;	#	14	K_gene_gntR
			1.0		;	#	15	𝛕_gene_gntR

			# translation kinetic limit parameters: P_venus
			1.0		;	#	16	K_P_venus
			1.0		;	#	17	𝛕_P_venus

			# translation kinetic limit parameters: P_gntR
			1.0		;	#	18	K_P_gntR
			1.0		;	#	19	𝛕_P_gntR

			# species degradation constants - 
			1.0		;	#	20	𝛳_mRNA_gntR
			1.0		;	#	21	𝛳_mRNA_venus
			1.0		;	#	22	𝛳_P_gntR
			1.0		;	#	23	𝛳_P_venus
		]


        # setup the parameter symbol index map -
        model_parameter_symbol_index_map = Dict{Symbol,Int}()
		model_parameter_symbol_index_map[:W_gene_venus] = 1
		model_parameter_symbol_index_map[:W_gene_venus_P_σ70] = 2
		model_parameter_symbol_index_map[:W_gene_venus_P_gntR] = 3
		model_parameter_symbol_index_map[:K_gene_venus_P_σ70] = 4
		model_parameter_symbol_index_map[:n_gene_venus_P_σ70] = 5
		model_parameter_symbol_index_map[:K_gene_venus_P_gntR] = 6
		model_parameter_symbol_index_map[:n_gene_venus_P_gntR] = 7
		model_parameter_symbol_index_map[:W_gene_gntR] = 8
		model_parameter_symbol_index_map[:W_gene_gntR_P_σ70] = 9
		model_parameter_symbol_index_map[:K_gene_gntR_P_σ70] = 10
		model_parameter_symbol_index_map[:n_gene_gntR_P_σ70] = 11
		model_parameter_symbol_index_map[:K_gene_venus] = 12
		model_parameter_symbol_index_map[:𝛕_gene_venus] = 13
		model_parameter_symbol_index_map[:K_gene_gntR] = 14
		model_parameter_symbol_index_map[:𝛕_gene_gntR] = 15
		model_parameter_symbol_index_map[:K_P_venus] = 16
		model_parameter_symbol_index_map[:𝛕_P_venus] = 17
		model_parameter_symbol_index_map[:K_P_gntR] = 18
		model_parameter_symbol_index_map[:𝛕_P_gntR] = 19
		model_parameter_symbol_index_map[:𝛳_mRNA_gntR] = 20
		model_parameter_symbol_index_map[:𝛳_mRNA_venus] = 21
		model_parameter_symbol_index_map[:𝛳_P_gntR] = 22
		model_parameter_symbol_index_map[:𝛳_P_venus] = 23


        # setup the inverse parameter symbol index map -
        inverse_model_parameter_symbol_index_map = Dict{Int, Symbol}()
		inverse_model_parameter_symbol_index_map[1] = :W_gene_venus
		inverse_model_parameter_symbol_index_map[2] = :W_gene_venus_P_σ70
		inverse_model_parameter_symbol_index_map[3] = :W_gene_venus_P_gntR
		inverse_model_parameter_symbol_index_map[4] = :K_gene_venus_P_σ70
		inverse_model_parameter_symbol_index_map[5] = :n_gene_venus_P_σ70
		inverse_model_parameter_symbol_index_map[6] = :K_gene_venus_P_gntR
		inverse_model_parameter_symbol_index_map[7] = :n_gene_venus_P_gntR
		inverse_model_parameter_symbol_index_map[8] = :W_gene_gntR
		inverse_model_parameter_symbol_index_map[9] = :W_gene_gntR_P_σ70
		inverse_model_parameter_symbol_index_map[10] = :K_gene_gntR_P_σ70
		inverse_model_parameter_symbol_index_map[11] = :n_gene_gntR_P_σ70
		inverse_model_parameter_symbol_index_map[12] = :K_gene_venus
		inverse_model_parameter_symbol_index_map[13] = :𝛕_gene_venus
		inverse_model_parameter_symbol_index_map[14] = :K_gene_gntR
		inverse_model_parameter_symbol_index_map[15] = :𝛕_gene_gntR
		inverse_model_parameter_symbol_index_map[16] = :K_P_venus
		inverse_model_parameter_symbol_index_map[17] = :𝛕_P_venus
		inverse_model_parameter_symbol_index_map[18] = :K_P_gntR
		inverse_model_parameter_symbol_index_map[19] = :𝛕_P_gntR
		inverse_model_parameter_symbol_index_map[20] = :𝛳_mRNA_gntR
		inverse_model_parameter_symbol_index_map[21] = :𝛳_mRNA_venus
		inverse_model_parameter_symbol_index_map[22] = :𝛳_P_gntR
		inverse_model_parameter_symbol_index_map[23] = :𝛳_P_venus


        # specific growth rate (default units: h^-1)
        μ = 0.0 # default units: h^-1

        # == DO NOT EDIT BELOW THIS LINE ========================================================== #
        problem_dictionary["initial_condition_array"] = initial_condition_array
        problem_dictionary["number_of_states"] = length(initial_condition_array)
        problem_dictionary["number_of_transcription_processes"] = number_of_transcription_processes
        problem_dictionary["number_of_translation_processes"] = number_of_translation_processes
        problem_dictionary["system_concentration_array"] = system_concentration_array
        problem_dictionary["biophysical_parameters_dictionary"] = biophysical_parameters_dictionary
        problem_dictionary["model_parameter_array"] = model_parameter_array
        problem_dictionary["model_parameter_symbol_index_map"] = model_parameter_symbol_index_map
        problem_dictionary["inverse_model_parameter_symbol_index_map"] = inverse_model_parameter_symbol_index_map
        problem_dictionary["stoichiometric_matrix"] = SM
        problem_dictionary["dilution_degradation_matrix"] = DM
        problem_dictionary["specific_growth_rate"] = μ

        # return -
        return problem_dictionary
        # ========================================================================================= #
    catch error
        throw(error)
    end
end
