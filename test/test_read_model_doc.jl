using JUGRN

# setup paths -
path_to_input_file = "/Users/jeffreyvarner/Desktop/julia_work/JUGRN.jl/test/data/Model.json"
path_to_output_dir = "/Users/jeffreyvarner/Desktop/julia_work/JUGRN.jl/test/generated_code"

# build the model object -
julia_model_object = build_julia_model_object(path_to_input_file, path_to_output_dir) |> check

# generate the ir dictionary -
ir_dictionary = read_model_document(julia_model_object) |> check
