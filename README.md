# General ODE systems integration library

This library is in development to provide a general API to call ODE
integration routines written in C. It is not meant to be a complete
database with all possible routines. The main goal is to provide the
most known and useful integrators with maximum performance possible
and a general API

## Application Building

To build the application, `cmake` is used. The following commands
build all application components:

```cmake
CC=gcc cmake -S . -B build
cmake --build build --target install
```

which must be executed from the `CMakeLists.txt` directory. In this
step, the `build`, `lib` and `bin` directories are created.

## Library usage

In directory `apps` some examples from references listed in the end
are implemented and the values rigorously compared with all decimal
places (when provided by the refs). In a first contact, reading the
executables in `apps` directory is highly recommended to understand
the workflow to use the library. Nevertheless, some of the assistant
routines are not employed in these executables to ease the lib usage
especially for large systems. The corresponding executables are
found in `bin` directory after running the `cmake` commands listed
above, with the same name of `apps` files. Before going over these
examples, the general workflow is explained in the following

### ODE system definition and derivatives evaluation

ODE systems must be taken in the standard form y' = f(x, y(x)) with
y' representing the derivative of y with respect to x. In this form
the most important element is the function implementation, which is
provided by the user. It must have the following signature:

- Real system `sys_der_func(RealODEInputParameters, double *)`
- Complex system `sys_der_func(ComplexODEInputParameters, double complex *)`

where `RealODEInputParameters` is a pointer to a `_RealODEInputParameters`
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
another struct, array of some data type and even another function.
This void pointer is also an explict parameter needed in the integrators,
and the user is reponsible to *cast* the pointer to a specific type
defined in the application client inside the function. An example
is provided in `apps/quinney_examples.c`

The app client do not need to use all `_RealODEInputParameters`
fields. The ODE system derivative signature was established in this
way since if all struct fields were explicitly passed as arguments,
some warnings might occur in compilation time if not all params
were used. Moreover, the choice of the pointer is to avoid multiple
struct copies and eventually improve performance.

The `_ComplexODEInputParameters` struct difference to the real case
is the substitution `Rarray -> Carray`, which correspond to `double *`
and `double complex *` pointers respectively.

Besides the function corresponding to y' = f(x, y(x)), the initial
condition must be provided by the application client as array with
the correponsing datatype of the ODE system (real or complex). For
multistep methods it is a bit more complicated to initiate and more
details are provided below.

### The workspace

After defining the function according to the signature above, before
going over the integration routines, a workspace must be provided in
order to avoid system calls for memory allocation during the ODE system
integration.

For instance, in Runge-Kutta methods, derivatives of the system are
evaluated in intermediate steps, and some workspaces arrays provide
a temporary place to store the output of function calls. In multistep
methods, in general, more than 1 step must be provided simultaneously
and a workspace to hold derivatives of several previous steps.

The following struct is used for runge-kutta methods

```
typedef struct{
    int system_size;
    Rarray
        work1,
        work2,
        work3,
        work4,
        work5;
} _RealWorkspaceRK;
```

and there is also `_ComplexWorkspaceRK` variant for complex systems.
Each one of these work arrays must have `system_size` elements. For
multistep method:

```
typedef struct{
    int
        ms_order,       /// number of previous steps required
        system_size;    /// number of equations in ODE system
    Rarray prev_der;    /// Hold all required previous derivatives
} _RealWorkspaceMS;
```

To provide all required previous derivatives the `prev_der` must
have size `ms_order * system_size`, which are concatenated from
the more recent to the older. More information about how to set
these previous steps in a single array are discussed below.
