## RecFMM build instructions with Meson

### Requirements

#### Build system

The updated recFMM port can be built using the [Meson][meson-home] build system
with the [Ninja][ninja-home] backend. These can be installed via `pip`:

    pip install meson
    pip install ninja

For more information on installing Meson and Ninja, refer to their respective
websites.

[meson-home]: https://mesonbuild.com/
[ninja-home]: https://ninja-org.build/

#### Compiler with Cilk support

The recFMM library supports shared-memory multi-threaded computations using the
Cilk C-language extensions and reducer hyperobjects.  To do that, you must use a
Cilk-supported C compiler.  RecFMM has been tested with the following C
compilers:

- [OpenCilk][opencilk-home] 1.0 (based on `clang` 10.0.1)
- GNU `gcc` 7.5.0
- Intel `icc` 19.1.2.254

We recommend using the newly released [OpenCilk][opencilk-home] compiler, as
Cilk support via the legacy Cilk Plus platform has been removed or deprecated in
other compilers.

A Fortran compiler is also required.  RecFMM has been tested with the
following Fortran compilers:

- GNU `gfortran` 7.5.0, 9.3.0
- Intel `ifort` 19.1.2.254

#### Memory allocation

Efficient memory allocation is supported via the Intel TBB scalable allocator
library, `tbbmalloc`, if that is available in your system.

### Building recFMM

#### Configuration

    CC=<c-compiler> FC=<fortran-compiler> meson path/to/build/directory

By default, Meson looks for dependencies using `pkg-config`.  To specify a
custom TBB installation directory, use

    CC=<c-compiler> FC=<fortran-compiler> meson -Dtbb_dir=path/to/tbb/directory path/to/build/directory

For additional configuration options, refer to the [Meson built-in
options][meson-builtin-options] web page.

[meson-builtin-options]: https://mesonbuild.com/Builtin-options.html

#### Compilation

    meson compile -C path/to/build/directory

#### Testing

    meson test -C path/to/build/directory

### Contributors

_Algorithm and library development:_
- Bo Zhang
- Jingfang Huang
- Nikos P. Pitsianis
- Xiaobai Sun

_Meson build config, OpenCilk support, bug fixes:_
- Alexandros-Stavros Iliopoulos
