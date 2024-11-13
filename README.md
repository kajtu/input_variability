# Input variability in a 4D flow MRI-based cardiovascular model
These are scripts to reproduce results in the manuscript *"Observer- and sequence variability in personalized 4D flow MRI-based cardiovascular models"* by Belén Casas Garcia, Kajsa Tunedal, Federica Viola, Gunnar Cedersund, Carl-Johan Carlhäll, Matts Karlsson, and Tino Ebbers.

The manuscript is available as a preprint at doi [10.1101/2024.06.13.597551](https://doi.org/10.1101/2024.06.13.597551)


If you use this implementation in your academic projects, please cite this paper.

# Input variability
The code for performing the analysis on inter-sequence and intra- and interobserver variability is provided in the scripts  `results_estimatedparameters.m ` and  `results_inputvars.m `.
Note that the original data is not provided here due to ethical restrictions.

# Sensitivity analysis
To perform the sensitivity analysis as described in the Supplementary, run `sensitivityanalysis.m `. 

# Requirements
The code was created with R2023a. Earlier Matlab versions might not be compatible with some of the scripts.

The model is implemented in the [AMICI toolbox](https://doi.org/10.1093/bioinformatics/btab227) in MATLAB for performing the sensitivity analysis.
To compile the model, MATLAB 2017b or earlier is needed, but to run the already compiled model any later matlab verison works. 
The provided compiled model is compiled on Windows, but will not work on Linux or macOS. 
To re-compile: run GenerateModels in the folder modelfiles. 
To compile the model, you need a valid C-compiler (such as xcode on Mac or MinGW on Windows. Run mex -setup to check if you have an installed compiler in matlab) and the MATLAB Symbolic Math Toolbox.


## Author
Kajsa Tunedal (kajsa.tunedal@liu.se) and Belén Casas Garcia.

## License
The MIT License (MIT)

Copyright (c) 2024 Kajsa Tunedal and Belén Casas Garcia 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

