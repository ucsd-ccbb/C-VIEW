#!/bin/bash

PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline

# Some string handling done here to remove spaces from result so that it
# can safely be passed as a variable to slurm's sbatch :-/
# tr ("translate") replaces all instances of first character w second character
# From https://linuxhint.com/newline_replace_sed/ :
# sed -z 's/\n/./g;s/.$/\n/' : The -z option is used to convert \n to the null
# character (\0). The content of the file is treated as a single line if it
# does not contain any null characters. The `sed` command will convert the
# newline into the null character and replace each \n with a comma by using
# the first search and replace pattern. Here, ‘g’ is used to globally search
# for \n. With the second search and replace pattern, the last comma will be
# replaced with \n.

cd $PIPELINEDIR && (git describe --tags && git log | head -n 1  && git checkout) | tr ' ' '_' | tr '\t' '_' | sed -z 's/\n/./g;s/.$/\n/'


