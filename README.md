# ENA_validation
Script to pre-screen gff3 for validation


# **Pre-requisites**

Python - 3.12 # Other versions may work but it has been tested with python 3.12

Python Libraries:
- collections
- argparse
- os


# **Installation**

This is a Python script and so once the script has been downloaded will run in any compatible python environment. **NB** If runnning on EI cluster the correct python environment will automatically be selected. If running outside of the institute then the shebang at the top of the python script (the first line) will have to be altered to reflect your python environment. A list of dependencies will follow shortly.

# How to run

If running on the EI cluster it is recommended you run it the following way:

interactive --mem 4G # depending on the size of the genome you are working with this should be sufficient however memory requirements may have to be increased for organisms with large genomes such as wheat.

remove_duplications_abutting.py -v -o {MyFile}.checked.gff3 {MyFile}.gff3 > {MyFile}.checked.log

-v = verbose – lots more info outputted to the log about what is being edited and removed
-o = what to call the output file

# **Licence**

This project is licensed under the GNU General Public License. See the LICENSE file for details.

# **Contact**

If you have any questions please create an issue or contact Camilla Ryan (camilla.ryan@earlham.ac.uk).
