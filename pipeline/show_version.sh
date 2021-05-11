#!/bin/bash

PIPELINEDIR=/shared/workspace/software/covid_sequencing_analysis_pipeline

cd $PIPELINEDIR && git describe --tags && git log | head -n 1 && git checkout

