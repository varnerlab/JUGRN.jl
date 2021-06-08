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


# calculate the u-variables -
function calculate_transcription_control_array(t::Float64, x::Array{Float64,1}, 
    parameter_dictionary::Dict{String,Any})::Array{Float64,1}

    # initialize -
    number_of_transcription_processes = parameter_dictionary["number_of_transcription_processes"]
    model_parameter_array = parameter_dictionary["model_parameter_array"]
    model_parameter_index_map = parameter_dictionary["model_parameter_symbol_index_map"]
    u_array = Array{Float64,1}()

    # local helper functions -
    f(x,K,n) = (x^n)/(K^n+x^n)
    u(A,R) = A/(1+A+R)

    # alias the state vector -
    gene_venus = x[1]
	mRNA_venus = x[2]
	P_venus = x[3]

    # alias the system species -
    system_array = parameter_dictionary["system_concentration_array"]
    RNAP = system_array[1]
	RIBOSOME = system_array[2]
	P_σ70 = system_array[3]


    # alias the system parameters -
    W_gene_venus = model_parameter_array[model_parameter_index_map[:W_gene_venus]]
	W_gene_venus_P_σ70 = model_parameter_array[model_parameter_index_map[:W_gene_venus_P_σ70]]
	K_gene_venus_P_σ70 = model_parameter_array[model_parameter_index_map[:K_gene_venus_P_σ70]]
	n_gene_venus_P_σ70 = model_parameter_array[model_parameter_index_map[:n_gene_venus_P_σ70]]
	K_gene_venus = model_parameter_array[model_parameter_index_map[:K_gene_venus]]
	𝛕_gene_venus = model_parameter_array[model_parameter_index_map[:𝛕_gene_venus]]
	K_P_venus = model_parameter_array[model_parameter_index_map[:K_P_venus]]
	𝛕_P_venus = model_parameter_array[model_parameter_index_map[:𝛕_P_venus]]
	𝛳_mRNA_venus = model_parameter_array[model_parameter_index_map[:𝛳_mRNA_venus]]
	𝛳_P_venus = model_parameter_array[model_parameter_index_map[:𝛳_P_venus]]


    # == CONTROL LOGIC BELOW ================================================================= #
    	
	# - gene_venus --------------------------------------------------------------------------------------------- 
	# gene_venus activation - 
	W = [
			W_gene_venus	;
			W_gene_venus_P_σ70	;
	]
	gene_venus_activator_array = Array{Float64,1}()
	push!(gene_venus_activator_array, 1.0)
	P_σ70_RNAP = RNAP*f(P_σ70,K_gene_venus_P_σ70,n_gene_venus_P_σ70)
	push!(gene_venus_activator_array, P_σ70_RNAP)
	A = sum(W.*gene_venus_activator_array)

	# gene_venus repression - 
	R = 0.0
	push!(u_array, u(A,R))


    # == CONTROL LOGIC ABOVE ================================================================= #

    # return -
    return u_array
end

# calculate the w-variables -
function calculate_translation_control_array(t::Float64, x::Array{Float64,1}, 
    parameter_dictionary::Dict{String,Any})::Array{Float64,1}
    
    # defualt: w = 1 for all translation processes -
    number_of_translation_processes = parameter_dictionary["number_of_translation_processes"]
    w = ones(number_of_translation_processes)
    
    # return -
    return w
end
