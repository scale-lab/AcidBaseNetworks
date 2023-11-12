Digital circuits and neural networks based on acid-base chemistry implemented by robotic fluid handling

Introduction:
This is the codebase and scripts of our paper: Digital circuits and neural networks based on acid-base chemistry implemented by robotic fluid handling (https://www.nature.com/articles/s41467-023-36206-8). 


1) Environment and Dependencies:
- Tested on Linux (Debian 9.3)
- Python 3.7
- Uses Python's virtualenv
- All needed Python packages are installable from dataset-gen/requirements.txt & simulations/requirements.txt as explained below.
- Pacakge dependencies are pretty light, hence, installation should be very fast depending on the bandwidth. 
- Every generation or simulation script takes around 1-10 mins depending on the resources and the datasize on average desktops.

The code is split into two parts:
- **dataset-gen/** contains the scripts that are used to download MNIST dataset, run the model training, binarization of the weights, and the dataset.
- **simulations/** contains the echo scripts and the simulation files to simulate the datasets generated from the previous part.

2) Dataset:
# Option A): Download the dataset:
To download the dataset files:
- Navigate to **simulations/**
- Run `python download_dataset.py`
- The script will download and extract the dataset to **simulations/simulation/datasets**

# Option B): Generate the dataset
It it recommended to use the dataset provided through option A, but if you prefer to generate the dataset by yourself:
- Navigate **to dataset-gen/**
- Create a Python virtual environment using: `virtualenv venv`
- Activate the environment: `source venv/bin/activate`
- Install the dependencies: `pip install -r requirements.txt`
- Finally, to train the model and generate the dataset files, run: `./generate_ds_imgs.sh`
- This will generate the needed data under **datasets/**

3) Running the simulations:
To run the simulations:
- Navigate to **simulations/**
- Create a Python virtual environment using: `virtualenv venv`
- Activate the environment: `source venv/bin/activate`
- Install the dependencies: `pip install -r requirements`
- Set the environment variable: `export PYTHONPATH=$(pwd):$(pwd)/chemcpupy`
- Set the variable **data_path** to the corresponding datasset path.
- Finally, to execute the simulations, run: `./sim_2_class_bin.sh`, `./sim_3_class_bin.sh`, or `sim_2_class_3bit.sh` to simulate the corresponding experiment from the paper.
- The scripts will generate accuracy files in the form of **acc\_\<bin or 3bit\>\_\<number of classe\>class\_\<image size\>.txt** that contains the accuracy metrics for the given simulations. The expected output should be similar to the results table in the manuscript.
- To clear the generated data, run `./clean.sh` to delete the generated files.

MNIST Dataset License:
Yann LeCun and Corinna Cortes hold the copyright of MNIST dataset, which is a derivative work from original NIST datasets. MNIST dataset is made available under the terms of the [Creative Commons Attribution-Share Alike 3.0 license](https://creativecommons.org/licenses/by-sa/3.0/).

License:
See "LICENSE.md" file

Citation:
If you find our work helpful in your research, please cite our paper:
@article{agiza2023digital,
  title={Digital circuits and neural networks based on acid-base chemistry implemented by robotic fluid handling},
  author={Agiza, Ahmed A and Oakley, Kady and Rosenstein, Jacob K and Rubenstein, Brenda M and Kim, Eunsuk and Riedel, Marc and Reda, Sherief},
  journal={Nature communications},
  volume={14},
  number={1},
  pages={496},
  year={2023},
  publisher={Nature Publishing Group UK London}
}
