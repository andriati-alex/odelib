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

with its fields explained by the comments. Note the void pointer
which provide a general API for the user freely tweak the system
as needed. This void pointer can take the address of anything, as
another struct, array of some data type and even another function.
This void pointer is also an explict parameter needed in the integrators,
and the user is reponsible to *cast* the pointer to a specific type
defined in the application client inside the function. An example
is provided in `apps/quinney_examples.c`

The app client does not need to use all `_RealODEInputParameters`
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
the corresponding datatype of the ODE system (real or complex). For
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

```c
typedef struct{
    unsigned int system_size;
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

```c
typedef struct{
    unsigned int
        ms_order,       /// number of previous steps required
        system_size;    /// number of equations in ODE system
    Rarray prev_der;    /// Hold all required previous derivatives
} _RealWorkspaceMS;
```

To provide all required previous derivatives, `prev_der` should have
at least size `ms_order * system_size`, with derivatives concatenated. 
Although, if the method is implicit one more chunk of size `system_size`
is required, due to implicit derivative evaluation of unknown step.
More information about how to set these previous steps in a single
array are discussed below.

### Assistant functions

After creating the structs described above, their fields must be set.
For that, some functions can help in this process, starting after the
fields have been set. These functions are:

- `alloc_real_rungekutta_wsarrays(_RealWorkspaceRK *)`
- `alloc_real_multistep_wsarray(_RealWorkspaceMS *)`

Which require the address of the struct created **WITH INTEGER FIELDS
ALREADY SET**. These functions will allocate all required workspace
memory before calling the integration routines without more labor.

Another way to set these workspace structs is using the automatic
allocation in HEAP with:

- `(_RealWorkspaceRK *) get_real_rungekutta_ws(int sys_size)`
- `(_RealWorkspaceMS *) get_real_multistep_ws(int sys_size, int ms_order)`

Where the returned pointer carry address of a struct ready to use.

### Short explanation about naming conventions

Just part of the functionality is been highlighted here, as how to
use for a real system of equation. However, there is also a complex
counterpart. In general, the following rules are used:

1. Camel case for structs datatypes and their pointers. Actual structs
   starts with an underscore and their pointers not.
2. Lower case with underscores for function names
3. In structs `Real` part of the (camel case) name can be changed to
   `Complex` for the equivalent data struct for a complex system
4. In function names, `real` can be changed to `cplx` for the equivalent
   complex case counterpart.
5. Also in function names, **in general**, there is a `rungekutta`
   which can be changed to `multistep` (does not hold for all)
6. Function names starting with `alloc` will allocate memory for
   arrays of real or complex numbers inside some struct. When
   starting with `free`, it means free all arrays in the struct
7. Function names starting with `get` will allocate structs in
   HEAP memory and return the pointer to it, with all internal
   fields set. Usually require integer arguments. When start
   with `destroy` do the opposite and free all internal pointers
   and the pointed struct itself.

### Main integration routines

With all the discussion of the previous sections, the integration
API is now easier to understand. First, to be clear in the following,
we can separate the problem in two parts, the known solution at
some grid point and the unknown solution at next grid point. We
can then use `y(x)` and `y(x+h)` or simply `y` and `ynext`.

To call the integrator to propagate one step, one needs to pass
the step size, the current grid point, the function reference to
evaluate derivatives, a pointer with extra arguments (as it is
void can be the address of anything), a workspace address (compatible
with the integrator routine), the current known step(s) required
and finally the array for output values with the result `ynext`

After calling the step integration the user is responsible to set
conditions for the next step, for instance, copying values from
`ynext -> y` and `x + h -> x`. It is also client's responsibility to
update extra arguments pointer if needed, as example for time
dependent parameters inside the system.

### A word about multistep integration

In multistep integration, the number of previous known steps required
depends on the multistep order defined in the workspace struct. Thus,
the current known step mentioned above `y` must actually be a collection
of `m` previous known steps, where `m` is the multistep order. Suppose we
want to obtain the `n+1` step solution, then we need `y_n ... y_n+1-m`.
Thus `y` must contain these steps concatenate, each one occupying a chunk
of `system_size` in the **VERY SAME ORDER PRESENTED HERE**, that is, from the
nearest known step to the oldest.

Moreover, the general multistep routine also require the set of coeff.
multiplying these steps, such that, it is possible to implement any
particular method, although the client is recommended to use the drivers
which implement particular schemes, such as Adams predictor correctors.

For those interested in the general multistep routine, it consumes two
extra parameters `a` and `b`, which are arrays of size `m + 1` with `m` the
multistep order (number of old steps used), and implies analytically in
the following equality solved numerically:

```
`a[0] * y_j+1 + a[1] * y_j + ... + a[m] * y_j+1-m =
			h * (b[0] * y'_j+1 + b[1] * y'_j + ... + b[m] * y'_j+1-m)`
```

where each one of these `yj` is a array (chunk of input known concat. steps `y`)

Finally, there is one more argument needed for multistep methods compared
to Runge-Kutta, whether the method is implicit or not. In case of implicit
method `ynext` itself become an input parameter, interpreted as the predictor
and the solution is obtained by iterating `niter` times, which is provided
as input parameter.

For more specific usage example, now its time to browse the `apps` files.



