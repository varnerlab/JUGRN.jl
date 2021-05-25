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
    gene_gntR = x[1]
	gene_venus = x[2]
	mRNA_gntR = x[3]
	mRNA_venus = x[4]
	P_gntR = x[5]
	P_venus = x[6]

    # alias the system species -
    system_array = problem_dictionary["system_concentration_array"]
    RNAP = system_array[1]
	RIBOSOME = system_array[2]
	σ70 = system_array[3]
	M_gluconate_c = system_array[4]


    # == CONTROL LOGIC BELOW ================================================================= #
    # venus_promoter: gene_venus -[RNAP]- mRNA_venus
	# venus_promoter activator set
	venus_promoter_activator_set = Array{Float64,1}()
	push!(venus_promoter_activator_set, 1.0)
	σ70_RNAP = RNAP*f(σ70,K,n)
	push!(venus_promoter_activator_set, σ70_RNAP)
	A = sum(W.*venus_promoter_activator_set)

	# venus_promoter repressor set
	venus_promoter_repressor_set = Array{Float64,1}()
	P_gntR_active = P_gntR*(1.0 - f(gluconate, K, n))
	push!(venus_promoter_repressor_set, P_gntR_active)
	R = sum(W.*venus_promoter_repressor_set)
	push!(u_array, u(A,R))
	# gntR_promoter: gene_gntR -[RNAP]- mRNA_gntR
	# gntR_promoter activator set
	gntR_promoter_activator_set = Array{Float64,1}()
	push!(gntR_promoter_activator_set, 1.0)
	σ70_RNAP = RNAP*f(σ70,K,n)
	push!(gntR_promoter_activator_set, σ70_RNAP)
	A = sum(W.*gntR_promoter_activator_set)
	R = 0
	push!(u_array, u(A,R))

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
