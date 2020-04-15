#!/bin/bash

conda env create -f environment.yml -p /home/user/anaconda3/envs/env_name
nglview enable
jupyter contrib nbextension install --user
jupyter-labextension install nglview-js-widgets

