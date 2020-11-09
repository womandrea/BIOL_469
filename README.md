# Any fans of phages?
Running phanotate, annotating what we get, and comparing it to those pesky bacterial annotaters.

## Dependencies 
Requires the pandas, BioPython, argParse, subprocess, os, and sys libraries to be accessible in the environment you run in. It would be best to create a separate ennvironment for this. You can do so by first ensuring you have miniconda install in your device, then running the following code:

```{bash}
conda create -n ENVIRONMENT_NAME
conda install pandas
conda install -c conda-forge/label/cf202003 argparse
conda install -c conda-forge biopython
```

All that is required is these libraries, the run_phanotate.py script, and for you to have Phanotate installed! You will need to update the file-path for Phanotate within the script. 
Good luck!!
