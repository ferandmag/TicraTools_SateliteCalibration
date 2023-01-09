ticra_analysis.py


can be used via command line by specifying the following arguments:

argvs[0] = ticra_analysis.py
argvs[1] = filename
argvs[2] = Te

filename (str)
Te (int)


parse and processes an output of the program ticra-tools-21.0.2.

creates a directory in the same place where the script is located called "Analysis YYYYY-MM-DD , HH-MM--SS , filename" in which the following files are stored: data (csv analysis data), figures (pdf with plots) and a series of csv matrices.

the file to be processed has to be located in the same folder where the script is located.


the file "CASS_planar_grid_B2_3_Casseg_MF" is left as an example.
