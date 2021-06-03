## Gene Regulatory Network Generator in Julia (JUGRN)

### Introduction ###
The ``JUGRN.jl`` package is a code generation system that transforms simple [JSON](https://www.json.org/json-en.html) descriptions of the connectivity of gene regulatory networks into model code written in the [Julia](http://julialang.org) programming language. ``JUGRN.jl`` has been used in the publications:

1. [Tasseff R, Jensen H, Congleton J, Dai W, Rogers K, Sagar A, Yen A and J. Varner (2017) An Effective Model of the Retinoic Acid Induced Differentiation Program, Sci Reports, 7:14327 doi:10.1038/s41598-017-14523-5](https://www.nature.com/articles/s41598-017-14523-5)
2. [Gould R, Bassen DM, Chakrabarti A, Varner JD and Butcher J (2016) Population Heterogeneity in the Epithelial to Mesenchymal Transition Is Controlled by NFAT and Phosphorylated Sp1. PLoS Comput Biol 12(12): e1005251. doi:10.1371/journal.pcbi.1005251](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005251)
3. [Adhikari A, Vilkhovoy M, Vadhin S, Lim HE and Varner JD (2020) Effective Biophysical Modeling of Cell Free Transcription and Translation Processes. Front. Bioeng. Biotechnol. 8:539081. doi: 10.3389/fbioe.2020.539081](https://www.frontiersin.org/articles/10.3389/fbioe.2020.539081/full)

to generate gene regulatory network code.

### Installation and Requirements
``JuGRN.jl`` is organized as a [Julia](http://julialang.org) package which 
can be installed in the ``package mode`` of Julia.

Start of the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/index.html) and enter the ``package mode`` using the ``]`` key (to get back press the ``backspace`` or ``^C`` keys). Then, at the prompt enter:

    (v1.1) pkg> add JUGRN

This will install the ``JUGRN.jl`` package and the other required packages. ``JUGRN.jl`` requires Julia 1.6.x and above.
``JUGRN.jl`` is open source. You can download this repository as a zip file, or clone or pull it by using the command (from the command-line):

	$ git pull https://github.com/varnerlab/JUGRN.jl.git

or

	$ git clone https://github.com/varnerlab/JUGRN.jl.git
