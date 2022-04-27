#!/bin/bash

PIPELINEDIR=/shared/workspace/software/cview

cd $PIPELINEDIR && git describe --tags && git log | head -n 1 && git checkout

