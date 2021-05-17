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
# THE SOFTWARE IS PROVIDED &quot;AS IS&quot;, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
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

    try

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
			1.0	;	#	σ70	units: µM
			1.0	;	#	M_gluconate_c	units: µM
		]

        # == DO NOT EDIT BELOW THIS LINE ======================================================= #
        problem_dictionary["initial_condition_array"] = initial_condition_array
        problem_dictionary["system_concentration_array"] = system_concentration_array
        return problem_dictionary
        # ====================================================================================== #
    catch error
        throw(error)
    end
end
