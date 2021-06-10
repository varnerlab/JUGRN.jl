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

function Balances(dx,x, parameter_dictionary,t)

    # system dimensions and structural matricies -
    number_of_states = parameter_dictionary["number_of_states"]
    DM = parameter_dictionary["dilution_degradation_matrix"]
    SM = parameter_dictionary["stoichiometric_matrix"]
    half_life_translation = parameter_dictionary["half_life_translation_capacity"]
     
    # calculate the TX and TL kinetic limit array -
    transcription_kinetic_limit_array = calculate_transcription_kinetic_limit_array(t,x,parameter_dictionary)
    translation_kinetic_limit_array = calculate_translation_kinetic_limit_array(t,x,parameter_dictionary)
    
    # calculate the TX and TL control array -
    u = calculate_transcription_control_array(t,x,parameter_dictionary)
    w = calculate_translation_control_array(t,x,parameter_dictionary)

    # calculate the rate of transcription and translation -
    rX = transcription_kinetic_limit_array.*u
    rL = translation_kinetic_limit_array.*w
    rV = [rX ; rL]

    # calculate the degradation and dilution rates -
    rd = calculate_dilution_degradation_array(t,x,parameter_dictionary)

    # compute the model equations -
    dxdt = SM*rV + DM*rd

    # package -
    for index = 1:number_of_states
        dx[index] = dxdt[index]
    end

    # extra species: global translation capacity -
	dx[4] = -(log(2)*(half_life_translation^-1))*x[4]
end
