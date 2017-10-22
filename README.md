# METHOD
-----------------------

## Multifluid Electromagneto-HydroDynamics 
---------------------------------------------

A three dimensional single- and multi-fluid EMHD solver, based in CUDA and C++. 

### Documentation
I have tried to maintain good documentation standards, but cant guarantee that everything you want will be here. If you are unsure of any of the functionality, find the respective header file for the class or function that you are curious about (this includes the abstract base classes for any derived classes). 

### Testing
We use the Google Test framework for unit testing. *Dont touch the GoogleTest directory!* Any tests are saved in the `Tests/Src` directory.

To build and run all tests, navigate to the `Tests` directory and run

  `make tests`
  
  
### Rootfinder
Some simulations will require the use of an N-dimensional footfinder, either for a (semi-) implicit time integrator or
for the conservative to primitive transformation. We have elected to use the CMINPACK library\*, and to use or implement any changes in the library, *cd* into the Cminpack directory and hit
  `make`
 to compile all the object files. Then, if the build was successful, for gods sake dont touch/look at this library again.


### Simulations
Simulations are run from the *main.cu* scripts. The Makefile in the Project folder links to the test directory, it is recommended that simulations are run with 

  `make all`
  
so that gtest will flag up any broken tests since the last change. Otherwise, simply use
  `make run`
to compile and run the simulation.


### Plotting Tools
The *Src* directory has a tool for interactively plotting the end state of a simulation. The `interactivePlot.py` script requires data to be saved after the simulation in the *Data*
folder. This is done automatically when using the SaveData class---call the class constructor with a pointer to the SimData class whose data you wish to save. E.g. enter 
  
  `SaveData save(&data);`
  
in *main* after the initial data has been evolved. Running the python script as main will load and store the data ready for plotting, and the easiest way to interact with the data is in a python environment such as spyder.



\* *found at https://github.com/devernay/cminpack, but due to the cryptic and poorly laid out package we have moved bits about and re-order various headers and includes. Most of the preprocessor stuff has been deleted (using cuda architechture will result in Cminpack reals defaulting to double precision), some functions have been excluded as they're not needed here, and now for any usage we just include the cminpack.h header file (as opposed to including the cuda scripts directly, which is horrid practice).*
