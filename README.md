# GSA Parallel Implementation

## Folder Structure

The main directory contains three folders:

- `src`: Contains the source code of the GSA algorithm.
- `PBS_scripts`: Contains scripts used to run the code on the cluster.
- `results`: Contains the results of tests performed on the cluster.

---

## Compiling the Code

To compile the code on the cluster, first load the required modules:

```bash
module load mpich-3.2
```

Inside the `src` folder, you will find the source code files and a `Makefile`.

To compile the code, navigate to the `src` directory and run:

```bash
make
```

This will create the executable file `gsa` in the `src` directory.

---

## Running the Code

From the `src` directory, the code can be executed using the following command:

```bash
mpiexec -n <number_of_processes> ./gsa <dimensions> <population_size> <iterations> <debug>
```

**Arguments:**

- `<number_of_processes>`: Number of processes to use in the parallel implementation.
- `<dimensions>`: Number of problem dimensions.
- `<population_size>`: Size of the population.
- `<iterations>`: Number of iterations to perform.
- `<debug>`: Debug mode:
  - `1`: Prints additional information (e.g., best solution found).
  - `0`: Prints results in CSV format.

---

## PBS Script

To run the code on the cluster, use the scripts located in the `PBS_scripts` folder.

**Important:** Update the path to the executable file in the script from:

```bash
/home/nicola.muraro/project_hpc/Parallel-implementation-of-Gravitational-Search-Algorithm/src/gsa
```

to the appropriate path on your system.

<br>
Then, submit the job to the queue with:

```bash
qsub <script_name>
```

The script will load the necessary modules and run the program with default arguments.

---

## Results

The `results` folder contains the outputs of the tests performed on the cluster.  
It includes a Python notebook used to average the results and generate tables and plots for analysis.
