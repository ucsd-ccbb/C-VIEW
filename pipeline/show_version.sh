#!/bin/bash

echo $(git describe --tags)
echo $(git log | head -n 1)
echo $(git checkout)

