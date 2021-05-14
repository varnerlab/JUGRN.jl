function build_julia_copyright_header_buffer(intermediate_dictionary::Dict{String,Any})
  
    # What is the current year?
    current_year = string(Dates.year(now()))
  
    # Get comment data from
    buffer = ""
    buffer *= "# ----------------------------------------------------------------------------------- #\n"
    buffer *= "# Copyright (c) $(current_year) Varnerlab\n"
    buffer *= "# Robert Frederick Smith School of Chemical and Biomolecular Engineering\n"
    buffer *= "# Cornell University, Ithaca NY 14850\n"
    buffer *= "#\n"
    buffer *= "# Permission is hereby granted, free of charge, to any person obtaining a copy\n"
    buffer *= "# of this software and associated documentation files (the \"Software\"), to deal\n"
    buffer *= "# in the Software without restriction, including without limitation the rights\n"
    buffer *= "# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
    buffer *= "# copies of the Software, and to permit persons to whom the Software is\n"
    buffer *= "# furnished to do so, subject to the following conditions:\n"
    buffer *= "#\n"
    buffer *= "# The above copyright notice and this permission notice shall be included in\n"
    buffer *= "# all copies or substantial portions of the Software.\n"
    buffer *= "#\n"
    buffer *= "# THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
    buffer *= "# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
    buffer *= "# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
    buffer *= "# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
    buffer *= "# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
    buffer *= "# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n"
    buffer *= "# THE SOFTWARE.\n"
    buffer *= "# ----------------------------------------------------------------------------------- #\n"
  
    # return -
    return buffer  
end