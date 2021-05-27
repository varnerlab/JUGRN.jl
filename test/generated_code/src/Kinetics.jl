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


function calculate_transcription_kinetic_limit_array(t::Float64, x::Array{Float64,1}, 
    problem_dictionary::Dict{String,Any})::Array{Float64,1}
    
    # initialize -
    kinetic_limit_array = Array{Float64,1}()
    system_array = problem_dictionary["system_concentration_array"]
    eX = problem_dictionary["biophysical_parameters_dictionary"]["transcription_elongation_rate"]       # default units: nt/s
    LX = problem_dictionary["biophysical_parameters_dictionary"]["characteristic_transcript_length"]    # default units: nt
    k_cat_characteristic = (eX/LX)

    # helper function -
    r(kcat, L_char, L, polymerase, 𝛕, K, species) = kcat*(L_char/L)*polymerase*(species/(𝛕*K+(1+𝛕)*species))

    # alias the model species -
    gene_gntR = x[1]
	gene_venus = x[2]
	mRNA_gntR = x[3]
	mRNA_venus = x[4]
	P_gntR = x[5]
	P_venus = x[6]

    # alias the system species -
    RNAP = system_array[1]
	RIBOSOME = system_array[2]
	P_σ70 = system_array[3]
	M_gluconate_c = system_array[4]


    # compute the transcription kinetic limit array -
    push!(kinetic_limit_array, r(k_cat_characteristic,LX,972,RNAP,𝛕_gene_venus,K_gene_venus,gene_venus))
	push!(kinetic_limit_array, r(k_cat_characteristic,LX,972,RNAP,𝛕_gene_gntR,K_gene_gntR,gene_gntR))

    # return -
    return kinetic_limit_array
end

function calculate_translation_kinetic_limit_array(t::Float64, x::Array{Float64,1}, 
    problem_dictionary::Dict{String,Any})::Array{Float64,1}

    # initialize -
    kinetic_limit_array = Array{Float64,1}()
    system_array = problem_dictionary["system_concentration_array"]
    eL = problem_dictionary["biophysical_parameters_dictionary"]["translation_elongation_rate"]         # default units: aa/s
    LL = problem_dictionary["biophysical_parameters_dictionary"]["characteristic_protein_length"]       # default units: aa
    k_cat_characteristic = (eL/LL)

    # helper function -
    r(kcat, L_char, L, polymerase, 𝛕, K, species) = kcat*(L_char/L)*polymerase*(species/(𝛕*K+(1+𝛕)*species))
    
    # alias the model species -
    gene_gntR = x[1]
	gene_venus = x[2]
	mRNA_gntR = x[3]
	mRNA_venus = x[4]
	P_gntR = x[5]
	P_venus = x[6]

    # alias the system species -
    RNAP = system_array[1]
	RIBOSOME = system_array[2]
	P_σ70 = system_array[3]
	M_gluconate_c = system_array[4]

    
    # compute the translation kinetic limit array -
    push!(kinetic_limit_array, r(k_cat_characteristic,LL,320,RIBOSOME,𝛕_P_venus,K_P_venus,P_venus))
	push!(kinetic_limit_array, r(k_cat_characteristic,LL,320,RIBOSOME,𝛕_P_gntR,K_P_gntR,P_gntR))

    # return -
    return kinetic_limit_array
end

function calculate_dilution_degradation_array(t::Float64, x::Array{Float64,1}, 
    problem_dictionary::Dict{String,Any})::Array{Float64,1}
    
    # initialize -
    degradation_dilution_array = Array{Float64,1}()
    
    # alias the model species -
    gene_gntR = x[1]
	gene_venus = x[2]
	mRNA_gntR = x[3]
	mRNA_venus = x[4]
	P_gntR = x[5]
	P_venus = x[6]

    # return -
    return degradation_dilution_array
end
