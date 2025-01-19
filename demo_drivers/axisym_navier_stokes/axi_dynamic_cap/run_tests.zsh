#!/usr/bin/zsh
find TEST -name parameters.dat | parallel --joblog tests.log --results log ./axi_dynamic_cap --parameters {}
