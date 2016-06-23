#!/bin/bash

dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

cd "$dir"

jupyter notebook ./gravitational_lensing.ipynb 


