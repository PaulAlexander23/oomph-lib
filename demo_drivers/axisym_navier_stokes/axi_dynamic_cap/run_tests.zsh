#!/usr/bin/zsh
find TEST -name parameters.dat | parallel --results log ./axi_dynamic_cap --parameters {}
