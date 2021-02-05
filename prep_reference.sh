#!/bin/bash

# prep reference
minimap2 -t THREADS -d NC045512.fas.mmi NC045512.fas
minimap2 -t THREADS -a -x map-ont REFERENCE_GENOME.FAS.MMI PRIMERS.FAS | samtools view -b -F 4 > PRIMERS.BAM
bedtools bamtobed -i PRIMERS.BAM > PRIMERS.BED