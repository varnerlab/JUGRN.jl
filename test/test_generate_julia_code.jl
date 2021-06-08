using JUGRNModelGenerator

# setup paths -
path_to_input_file = "/Users/jeffreyvarner/Desktop/julia_work/JUGRNModelGenerator.jl/test/data/Venus.json"
path_to_output_dir = "/Users/jeffreyvarner/Desktop/julia_work/JUGRNModelGenerator.jl/test/generated_code"

# build the model object -
julia_model_object = build_julia_model_object(path_to_input_file, path_to_output_dir) |> check

# try to generate code -
result = make_julia_model(julia_model_object)

