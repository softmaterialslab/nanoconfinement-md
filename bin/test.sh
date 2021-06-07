#!/bin/bash

python3 -m venv ENV

source ENV/bin/activate

pip install hublib
pip install jupyterlab
pip install testbook
pip install widgetsnbextension
pip install svgwrite
pip install pypng
pip install pytest
pip install scipy
pip install tensorflow
pip install keras
pip install pandas
pip install scikit-learn
pip install matplotlib
pip install ipywidgets
pip install pandas
jupyter nbextension enable --py widgetsnbextension

export ENVIRON_CONFIG_DIRS=","
export PATH=$PATH:$PWD
pytest test || exit 1
