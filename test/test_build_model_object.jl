using JUGRN

# setup paths -
path_to_input_file = "./test/data/Model.json"
path_to_output_dir = "./test/generated_code"

# build the model object -
julia_model_object = build_julia_model_object(path_to_input_file, path_to_output_dir)