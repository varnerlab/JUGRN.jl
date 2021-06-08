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

# include -
include("Include.jl")

function solve_dynamic_problem(time_start::Float64, time_stop::Float64, time_step::Float64, 
    parameter_dictionary::Dict{String,Any}, varargs...)

    # Get required stuff from the problem struct -
    time_span = (time_start,time_stop)
    initial_condition_array = parameter_dictionary["initial_condition_array"];

    # build problem object -
    problem_object = ODEProblem(Balances, initial_condition_array, time_span, parameter_dictionary)

    # solve -
    solution = solve(problem_object, AutoTsit5(Rosenbrock23(autodiff=false)), reltol=1e-8,abstol=1e-8)

    # pull solution apart -
    T = solution.t

    # initialize the state array -
    number_of_times_steps = length(T)
    number_of_states = length(initial_condition_array)
    X = zeros(number_of_times_steps,number_of_states)
    for step_index=1:number_of_times_steps

        # grab the solution 
        soln_array = solution.u[step_index]
        for state_index = 1:number_of_states
            X[step_index, state_index] = soln_array[state_index]
        end
    end

    # return -
    return (T,X)
end

function solve_static_problem(parameter_dictionary::Dict{String,Any}, varargs...)

    # Get the initial conditions -
    initial_condition_array = problem["initial_condition_array"];

    # calculate the steady-state soln -
    steady_state_prob = SteadyStateProblem(Balances, initial_condition_array, parameter_dictionary)
    steady_state_soln  = solve(steady_state_prob, DynamicSS(AutoTsit5(Rosenbrock23(autodiff=false)); abstol=1e-8, reltol=1e-6, tspan=Inf))

    # return -
    return abs.(steady_state_soln.u)
end
