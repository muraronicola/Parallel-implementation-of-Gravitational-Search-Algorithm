# GSA parallel implementation


## Folder structure
In the principal folder, there are tree folders:
- `src`: contains the source code of the GSA algorithm.
- `PBS_scripts`: contains the scripts used to run the code on the cluster.
- `results`: contains the results of the tests performed on the cluster.

<br>


## Compiling the code
For compiling the code on the cluster, first we need to load the required modules:
```bash
module load mpich-3.2
```

Inside the `src` folder, there are the source code files and a Makefile.
The Makefile can be used to compile the code, by simply entering the directory `src` and running the command:
```bash
make
```
This will create the executable file `gsa` in the `src` directory.

<br>


## Running the code
To run the code, from the `src` directory, we can use the following command:
```bash
mpiexec -n <number_of_processes> ./gsa <dimensions> <population_size> <iterations> <debug>
```
Where:
- `<number_of_processes>`: Number of processes to be used in the parallel implementation.
- `<dimensions>`: Number of dimensions of the problem.
- `<population_size>`: Size of the population.
- `<iterations>`: Number of iterations to be performed.
- `<debug>`: Debug mode. If set to 1, the code will print additional information, like the best solution at the end. Otherwise, if set to 0, the code will print the results in a .csv style file.

<br>


## PBS script
In order to run the code on the cluster, we can use the scripts present in the folder `PBS_scripts`.
Is necessary to change the path of the executable file from `/home/nicola.muraro/project_hpc/GSA/src/gsa` to the path of the executable file in your system.
Then, by simply running the command:
```bash
qsub <script_name>
```
we can submit the job to the queue. The script will take care of loading the required modules and running the programm with the default arguments.

<br>

## Results
The folder `results` contains the results of the tests performed on the cluster. There is a python notebook that has been used to average the results and obtain the tables and plots.
