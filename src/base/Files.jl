"""
    read_model_document(path_to_file::String;
            strip_comments::Bool = true)::Array{String,1}

Read each line of the Network.vff input file and extract contents of the file. Comments in the file (beginning by `//` symbol) are excluded by default, which can be overridden by the user.

Input arguments:
`path_to_file::String` - user-specified path where the Network.vff input file can be found and parsed
`strip_comments::Bool` - if true, comment lines (beginning by ‘//’ symbol) are excluded by the parser (optional).

Output arguments:
`buffer::Array{String,1}` - single array holding information extracted from the Network.vff file.

"""
function read_model_document(path_to_file::String;
    strip_comments::Bool = true)::Array{String,1}

    # initialize -
    buffer = String[]

    # Read in the file -
    open("$(path_to_file)", "r") do file
        for line in eachline(file)

            # exclude comments -
            if (occursin("//",line) == false && isempty(line) == false && strip_comments == true)
                +(buffer,line)
            end
        end
    end

    # return -
    return buffer
end

function transfer_distribution_file(path_to_distribution_files::String,
                                      input_file_name_with_extension::String,
                                      path_to_output_files::String,
                                      output_file_name_with_extension::String)

    # Load the specific file -
    # create src_buffer -
    src_buffer::Array{String} = String[]

    # check - do we have the file path?
    if (isdir(path_to_output_files) == false)
        mkpath(path_to_output_files)
    end

    # path to distrubtion -
    path_to_src_file = path_to_distribution_files*"/"*input_file_name_with_extension
    open(path_to_src_file,"r") do src_file
        for line in eachline(src_file)

            # need to add a new line for some reason in Julia 0.6
            new_line_with_line_ending = line*"\n"
            push!(src_buffer,new_line_with_line_ending)
        end
    end

    # Write the file to the output -
    path_to_program_file = path_to_output_files*"/"*output_file_name_with_extension
    open(path_to_program_file, "w") do f
        for line in src_buffer
            write(f,line)
        end
    end
end

function transfer_distribution_files(path_to_distribution_files::String,
                                      path_to_output_files::String,
                                      file_extension::String)


    # Search the directory for src files -
    # load the files -
    searchdir(path,key) = filter(x->contains(x,key),readdir(path))

    # build src file list -
    list_of_src_files = searchdir(path_to_distribution_files,file_extension)

    # check - do we have the file path?
    if (isdir(path_to_output_files) == false)
        mkpath(path_to_output_files)
    end

    # go thru the src file list, and copy the files to the output path -
    for src_file in list_of_src_files

        # create src_buffer -
        src_buffer::Array{String,1} = String[]

        # path to distrubtion -
        path_to_src_file = path_to_distribution_files*"/"*src_file
        open(path_to_src_file,"r") do src_file
            for line in eachline(src_file)

                # need to add a new line for some reason in Julia 0.6
                new_line_with_line_ending = line*"\n"
                push!(src_buffer,new_line_with_line_ending)
            end
        end

        # Write the file to the output -
        path_to_program_file = path_to_output_files*"/"*src_file
        open(path_to_program_file, "w") do f
            for line in src_buffer
                write(f,line)
            end
        end
    end
end

function write_program_components_to_disk(file_path::String, set_of_program_components::Set{NamedTuple})

    # check - do we have teh file path?
    if (isdir(file_path) == false)
      mkpath(file_path)
    end

    # go through each component, and dump the buffer to disk -
    for program_component in set_of_program_components

      # We switch on type -
      filename = program_component.filename
      component_type = program_component.component_type
      if (component_type == :buffer)

          # get the data -
          program_buffer = program_component.buffer

          # build the path -
          path_to_program_file = file_path*"/"*filename

          # Write the file -
          outfile = open(path_to_program_file, "w")
          write(outfile,program_buffer);
          close(outfile);

      elseif (component_type == :matrix || component_type == :vector)

          # get the matrix -
          program_matrix = program_component.matrix

          # build the path -
          path_to_program_file = file_path*"/"*filename

          # write the file -
          writedlm(path_to_program_file, program_matrix)

      else
          error("unsupported program component type: $(component_type)")
      end
    end
end

function move_existing_project_at_path(path_to_existing_project::String)::Bool

    # we are getting called *if* we already know there is a dir conflict -
    # if this is getting called, we have an existing dir where the user wants to write code.
    # we need then create a new dir called *.0, and mv the offending dir to this location?
    # return true if this worked, otherwise false -

    # parent and child dir that we are generating into -
    parent_dir = dirname(path_to_existing_project)
    child_dir = basename(path_to_existing_project)
    destination_path = ""

    # current backup index  -
    current_backup_index = 0

    # do we already have the destination?
    loop_flag = true
    while loop_flag

         # make a destination path -
        destination_path = joinpath(parent_dir,"$(child_dir).$(current_backup_index)")

        # we don't have this dir, we are done -
        if (isdir(destination_path) == false)
            loop_flag = false
        end

        # ok, looks like we already have this dir, update the counter -
        current_backup_index = current_backup_index + 1
    end

    # mv -
    mv(path_to_existing_project, destination_path)

    # check -
    if (isdir(destination_path) == false)
        return false
    end

    return true
end


"""
    generate_default_project(path_to_project_dir::String)::VLResult

Generates a default project structure which contains an empty model file and Defaults.toml file.

Input arguments:
`path_to_project_file::String` - path to where you want model code to be generated

Output arguments:
None
"""
function generate_default_project_file(path_to_project_file::String)::VLResult

    # ok, if we get here, then we have a clean place to generate the default project structure -
    # We need to two things for a project, the defaults file, and a blank network file with all the sections -
    filename = basename(path_to_project_file)
    path_to_output_files = dirname(path_to_project_file)

    # Transfer distrubtion files to the output -> these files are shared between model types -
    # TODO: replace w/joinpath?
    transfer_distribution_file("$(_PATH_TO_SRC)/distribution/dynamic", filename, path_to_output_files, filename)

    # ok, so return nothing if everything worked ...
    return VLResult(nothing)
end