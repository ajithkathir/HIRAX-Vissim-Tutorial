# HIRAX-Vissim-Tutorial
This is a repository for HIRAX visibility simulation tutorial meant for the usage of HIRAX collaborators only. 

# Instructions to setup

Make sure you have python3.10. Optionally, you can create a python3.10 environment using conda if you are dealing with multiple python versions. 

After cloning the repo, use the install.sh file to install necessary packages by running 
```bash
chmod +x install.sh
./install.sh 
```
This will ensure that you have a virtual environment where all packages are installed. 

In the terminal run :

```bash
source hiraxvenv/bin/activate
```
to activate the created virtual environment. 


Once the packages are installed and venv activated, we will use Jupyterlab or Notebook for the tutorial session. 
I have provided a sample skyimage file (h5) for this tutorial shared in googledrive folder.
Link: https://drive.google.com/drive/folders/1HTaJLolP4IzMQ6egkChvSUzNj93bqJnl?usp=drive_link

Download and place it in the tutorial folder. 

Also create a directory called /.casa/data in your home directory. For example:

```bash
mkdir /home/ajith/.casa/data
```
so that the CASA downloads latest measures data from its ftp servers, when it is run for the first time. 
