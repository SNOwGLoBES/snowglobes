# SNOwGLoBES
SuperNova Observatories with GLoBES

## Dependencies

Make sure you have **gcc** and **GNU gls** installed (GLoBES install prerequisite) <br>
  * **ROOT** must be installed ! 
  * **gcc** can be installed through the ROOT prereqs 
  * **GNU gls** can be intalled through the ROOT Optional prereqs
     *  See: https://root.cern/install/dependencies/ 

## Installation GLoBES
This software requires the GLoBES libraries to be installed.<br>
_See: https://www.mpi-hd.mpg.de/personalhomes/globes/_<br>

Make a parent dir for GLoBES/SNOwGLoBES: 
```bash
mkdir SN_stuff 
``` 

Move into SN_stff: 
```bash
cd SN_stff 
``` 

Download GLoBES: 
```bash
wget https://www.mpi-hd.mpg.de/personalhomes/globes/download/globes-3.2.18.tar.gz 
```

Untar the file: 
```bash 
tar xvfz globes-3.2.18.tar.gz
```

Go into the GLoBES dir and configure the install :
```bash
cd globes-3.2.18
./configure 
``` 

Run make and make install: 
```bash
make
make install
```

run ldconfig:
```bash
sudo ldconfig
```

## Installation SNOwGLoBES
Clone this repository:
```bash
git clone https://github.com/SNOwGLoBES/snowglobes.git
```

Set environment variables `$SNOWGLOBES` and `$GLB_DIR` to the directories where these are or will be installed
```bash
export GLB_DIR=/home/[user]/SN_stuff/globes-3.2.18
export SNOWGLOBES=/home/[user]/SN_stuff/snowglobes
```

Change into `$SNOWGLOBES/src` and type
```bash
make
make install
```

## Usage
```bash
./supernova.pl <run_mode> <flux_name> <channel_name> <detector>

# for example:
./supernova.pl 0 livermore lead halo1
```
