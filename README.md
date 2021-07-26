# General ODE systems integration library

This library is in development to provide a general API to call ODE
integration routines written in C. It is not meant to be a complete
database with all possible routines. The main goal is to provide the
most known and useful integrators with maximum performance possible
and a general API

## Application Building

To build the application the `cmake` tool is required. The following
commands are used to build all application components:

```cmake
CC=gcc cmake -S . -B build
cmake --build build --target install
```

which must be ran from the same `CMakeLists.txt` directory.
In this step, the `build`, `lib` and `bin` directories are created,
from which the important file is the `libode` in `lib` directory.

## Library usage

In directory `apps` some examples from references listed in the end
are implemented and the values rigorously compared with all decimal
places (when provided by the refs). In a first contact, reading the
executables in `apps` directory is highly recommended to understand
the workflow to use the library. The corresponding executables are
found in `bin` directory with the same name of `apps` files. Before
going over the example files the workflow is explained in the following

### ODE system definition and derivatives evaluation

ODE systems must be taken in the standard form y' = f(x, y(x)) with
y' representing the derivative of y with respect to x. In this form
the most important element is the function implementation, which is
provided by the user, must have the following signature:

- Real system `sys_der_func(RealODEInputParameters, double *)`
- Complex system `sys_der_func(ComplexODEInputParameters, double complex *)`

where, `RealODEInputParameters` is a pointer to a `_RealODEInputParameters`
struct which is defined as:

```c
/** Struct with input parameters for derivatives computation */
typedef struct{
    unsigned int system_size;   /// number of equations in the system
    double x;                   /// grid point of the known solution
    Rarray y;                   /// function values at `x`
    void * extra_args;          /// user-defined external arguments
} _RealODEInputParameters;

typedef _RealODEInputParameters * RealODEInputParameters;
```

with its fields exaplained by the comments. Note the void pointer
which provide a general API for the user freely tweak the system
as needed. This void pointer can take the address of anything, as
another struct for instance or simple array. This void pointer is
also an explict parameter needed in the integrators, and the user
is reponsible to *cast* the pointer to a specific type defined in
the application client(see `apps/quinney_examples.c` for example)

The app client is not required to use all `_RealODEInputParameters`
fields. The ODE system derivative signature was established in this
way since if each of these struct fields were explicitly passed as
arguments, some warning occur in compilation time if not all params
are used. Moreover, the choice of the pointer is to avoid multiple
struct copies and improve performance.

The `_ComplexODEInputParameters` only diff by `Rarray -> Carray`
which correspond to `double *` and `double complex *` pointers.

Besides the function corresponding to y' = f(x, y(x)), the initial
condition must be provided by the application client as array with
the correponsing datatype of the ODE system (real or complex)

### The workspace

After defining the function according to the signature above, before
going over the step-integrator routines, a workspace must be provided.
The workspace is introduce to avoid system calls for memory allocation
during the ODE system integration.
