# EvoDevo

Working draft of a model whereby the fitness of individuals from an evolving population is determined by their interactions when grouped to make small random Boolean networks

This code uses the json and h5 libraries via the morphologica library and you will need to install any dependencies for morphologica by following the instructions at https://github.com/ABRG-Models/morphologica.git.

Once you've done that clone this repo, and in this directory also clone morphologica:

```
git clone https://github.com/ABRG-Models/morphologica.git
```

Then you can build the project: 

```
mkdir build
cd build
cmake ..
make
cd ..
```

To run the simulation:

```
./build/evodevo config.json logs
```

(optionally append a seed for the random number generator)

When the simulation finishes, and logs/out.h5 has been created, you can do:

```
python analysis.py
```

(optionally provide the path to an alternative version of logs/out.h5)






