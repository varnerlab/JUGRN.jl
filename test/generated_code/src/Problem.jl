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
    system_type_flag = 	:CF_PURE

    try

        # load the stoichiometric_matrix (SM) and degradation_dilution_matrix (DM) -
        SM = readdlm("./network/Network.dat")
        DM = readdlm("./network/Degradation.dat")

        # build the species initial condition array -
        initial_condition_array = [
			5.0	;	#	1	gene_gntR	units: nM
			5.0	;	#	2	gene_venus	units: nM
			0.0	;	#	3	mRNA_gntR	units: nM
			0.0	;	#	4	mRNA_venus	units: nM
			0.0	;	#	5	P_gntR	units: nM
			0.0	;	#	6	P_venus	units: nM
		]

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
			0.001	;	#	1	W_gene_venus
			1.0		;	#	2	W_gene_venus_P_σ70
			1.0		;	#	3	W_gene_venus_P_gntR

			0.001	;	#	4	W_gene_gntR
			1.0		;	#	5	W_gene_gntR_P_σ70

			1.0		;	#	6	K_gene_venus
			1.0		;	#	7	K_gene_gntR

			1.0		;	#	8	𝛕_gene_venus
			1.0		;	#	9	𝛕_gene_gntR

			1.0		;	#	10	K_P_venus
			1.0		;	#	11	K_P_gntR

			1.0		;	#	12	𝛕_P_venus
			1.0		;	#	13	𝛕_P_gntR
		]


        # setup the parameter symbol - index map -
        

        # == DO NOT EDIT BELOW THIS LINE ========================================================== #
        problem_dictionary["initial_condition_array"] = initial_condition_array
        problem_dictionary["system_concentration_array"] = system_concentration_array
        problem_dictionary["biophysical_parameters_dictionary"] = biophysical_parameters_dictionary
        problem_dictionary["model_parameter_array"] = model_parameter_array
        problem_dictionary["model_parameter_symbol_index_map"] = model_parameter_symbol_index_map
        problem_dictionary["stoichiometric_matrix"] = SM
        problem_dictionary["dilution_degradation_matrix"] = DM

        # return -
        return problem_dictionary
        # ========================================================================================= #
    catch error
        throw(error)
    end
end
