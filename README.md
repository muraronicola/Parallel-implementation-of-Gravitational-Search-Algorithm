# GSA parallel implementation

## Compiling the code
In the principal folder, we can see all the files used in the implementation of the GSA algorithm.

In the case of compliling the code on the cluster, first we need to load the required modules:
```bash
module load mpich-3.2
```

There is a Makefile that can be used to compile the code, by simply running the command:
```bash
make
```
This will create the executable file `gsa` in the current directory.


## Running the code
To run the code, we can use the following command:
```bash
mpiexec -n <number_of_processes> ./gsa <dimentions> <population_size> <iterations> <debug>
```
Where:
- `<number_of_processes>`: Number of processes to be used in the parallel implementation.
- `<dimentions>`: Number of dimensions of the problem.
- `<population_size>`: Size of the population.
- `<iterations>`: Number of iterations to be performed.
- `<debug>`: Debug mode. If set to 1, the code will print additional information, like the best solution at the end.


## PBS script
In order to run the code on the cluster, we can use the scripts present in the folder `scripts`.
By simply running the command:
```bash
qsub <script_name>
```
We can submit the job to the queue. The script will take care of loading the required modules and running the programm with the default arguments.


## Results
The folder `results` contains the results of the tests performed on the cluster. There are two python notebook that has been used to average the results and obtain the tables and plots.
