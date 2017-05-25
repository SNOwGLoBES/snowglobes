# SNOwGLoBES
SuperNova Observatories with GLoBES

## Dependencies
This software requires the GLoBES libraries to be installed.<br>
_See: https://www.mpi-hd.mpg.de/personalhomes/globes/_

## Installation
Clone this repository:
```bash
git clone https://github.com/SNOwGLoBES/snowglobes.git
```

Set environment variables `$SNOWGLOBES` and `$GLB_DIR` to the directories where these are or will be installed
```bash
export GLB_DIR=/path/to/GLoBES/
export SNOWGLOBES=/path/to/SNOwGLoBES
```

Change into `$SNOWGLOBES/src` and type
```bash
make
make install
```

## Usage
```bash
./supernova.pl <flux_name> <channel_name> <detector>

# for example:
./supernova.pl livermore lead halo1
```
