                         Intel(R) Fortran Compiler Help
                         ==============================

  usage: ifort [options] file1 [file2 ...]

     where options represents zero or more compiler options

     fileN is a Fortran source (.f .for .ftn .f90 .fpp .F .FOR .F90 .i .i90),
     assembly (.s .S), object (.o), static library (.a), or other 
     linkable file

     Commonly used options may be placed in the ifort .cfg file.

   Some options listed are only available on a specific system
   i32    indicates the feature is available on systems based on IA-32
          architecture
   i64em  indicates the feature is available on systems using Intel(R) 64
          architecture

                             Compiler Option List
                             --------------------

Optimization
------------

-O1       optimize for maximum speed, but disable some optimizations which
          increase code size for a small speed benefit
-O2       optimize for maximum speed (DEFAULT)
-O3       optimize for maximum speed and enable more aggressive optimizations
          that may not improve performance on some programs
-O        same as -O2
-Os       enable speed optimizations, but disable some optimizations which
          increase code size for small speed benefit 
-O0       disable optimizations
-fast     enable -xHOST -O3 -ipo -no-prec-div -static
          options set by -fast cannot be overridden with the exception of
          -xHOST, list options separately to change behavior
-fno-alias
          assume no aliasing in program
-fno-fnalias
          assume no aliasing within functions, but assume aliasing across calls
-nolib-inline
          disable inline expansion of intrinsic functions

Code Generation
---------------

-x<code1>[,<code2>,...]
          generate specialized code to run exclusively on processors
          indicated by <code> as described below
            Host generate instructions for the highest instruction set and
                 processor available on the compilation host machine
            SSE2 Intel Pentium 4 and compatible Intel processors.  Enables new
                 optimizations in addition to Intel processor-specific
                 optimizations
            SSE3    Intel(R) Core(TM) processor family with Streaming SIMD
                    Extensions 3 (Intel(R) SSE3) instruction support
            SSE3_ATOM Can generate MOVBE instructions for Intel processors and
                    can optimize for the Intel(R) Atom(TM) processor.
            SSSE3   Intel(R) Core(TM)2 processor family with Supplemental
                    Streaming SIMD Extensions 3 (SSSE3)
            SSE4.1  Intel(R) 45nm Hi-k next generation Intel Core(TM)
                    microarchitecture with support for SSE4 Vectorizing
                    Compiler and Media Accelerator instructions
            SSE4.2  Can generate Intel(R) SSE4 Efficient Accelerated String
                    and Text Processing instructions supported by Intel(R)
                    Core(TM) i7 processors. Can generate Intel(R) SSE4 
                    Vectorizing Compiler and Media Accelerator, Intel(R) SSSE3,
                    SSE3, SSE2, and SSE instructions and it can optimize for
                    the Intel(R) Core(TM) processor family.
            AVX     Enable Intel(R) Advanced Vector Extensions instructions
-ax<code1>[,<code2>,...]
          generate code specialized for processors specified by <codes>
          while also generating generic IA-32 instructions.  
          <codes> includes one or more of the following:
            SSE2 Intel Pentium 4 and compatible Intel processors.  Enables new
                 optimizations in addition to Intel processor-specific
                 optimizations
            SSE3    Intel(R) Core(TM) processor family with Streaming SIMD
                    Extensions 3 (Intel(R) SSE3) instruction support
            SSSE3   Intel(R) Core(TM)2 processor family with Supplemental
                    Streaming SIMD Extensions 3 (SSSE3)
            SSE4.1  Intel(R) 45nm Hi-k next generation Intel Core(TM)
                    microarchitecture with support for Streaming SIMD
                    Extensions 4 (Intel(R) SSE4) Vectorizing
                    Compiler and Media Accelerator instructions
            SSE4.2  Can generate Intel(R) SSE4 Efficient Accelerated String
                    and Text Processing instructions supported by Intel(R)
                    Core(TM) i7 processors. Can generate Intel(R) SSE4 
                    Vectorizing Compiler and Media Accelerator, Intel(R) SSSE3,
                    SSE3, SSE2, and SSE instructions and it can optimize for
                    the Intel(R) Core(TM) processor family.
            AVX     Enable Intel(R) Advanced Vector Extensions instructions
-tune <keyword>
          pn1 - optimize for Pentium(R) processor (i32 only)
          pn2 - optimize for Pentium(R) Pro, Pentium(R) II, and Pentium(R) III
                processors (i32 only)
          pn3 - same as pn2 (i32 only)
          pn4 - optimize for Pentium(R) 4 processor (DEFAULT)
-arch <keyword>
          pn1    - optimize for Pentium(R) processor (i32 only)
          pn2    - optimize for Pentium(R) Pro, Pentium(R) II, and 
                   Pentium(R) III processors (i32 only)
          pn3    - same as pn2 (i32 only)
          pn4    - optimize for Pentium(R) 4 processor (DEFAULT)
          SSE    - Intel Pentium III and compatible Intel processors (i32 only)
          SSE2   - Intel Pentium 4 and compatible Intel processors.  Enables
                   new optimizations in addition to Intel processor-specific
                   optimizations (i32 only)
          SSE3   - Intel(R) Core(TM) processor family.  Code is expected to run
                   properly on any processor that supports SSE3, SSE2 and SSE
                   instruction sets
          SSSE3  - Intel(R) Core(TM)2 processor family with Supplemental 
                   Streaming SIMD Extensions 3 (SSSE3)
          SSE4.1 - Intel(R) 45nm Hi-k next generation Intel Core(TM)
                   microarchitecture with support for Streaming SIMD
                   Extensions 4 (Intel(R) SSE4) Vectorizing Compiler
                   and Media Accelerator instructions
-mcpu=<cpu>
          same as -mtune=<cpu>
-mtune=<cpu>
          optimize for a specific <cpu>
            pentium3  - optimize for Pentium(R) III processors
            pentium4  - optimize for Pentium(R) 4 processor (DEFAULT)
-march=<cpu>
          generate code exclusively for a given <cpu>
            pentium3  - streaming SIMD extensions
            pentium4  - Pentium(R) 4 New Instructions
-msse3    generate code for Intel(R) Core(TM) Duo processors, Intel(R) Core(TM)
          Solo processors, Intel Pentium 4 and compatible Intel processors with
          Streaming SIMD Extensions 3 (Intel(R) SSE3) instruction support
-mssse3   Intel(R) Core(TM)2 processor family with Supplemental Streaming 
          SIMD Extensions 3 (SSSE3)
-msse4.1  Intel(R) 45nm Hi-k next generation Intel Core(TM) microarchitecture
          with support for Streaming SIMD Extensions 4 (Intel(R) SSE4) 
          Vectorizing Compiler and Media Accelerator instructions
-minstruction=<keyword>
          Refine instruction set output for the selected target processor

            [no]movbe  - Do/do not generate MOVBE instructions with SSE3_ATOM
                          (requires -xSSE3_ATOM)
           
-f[no-]omit-frame-pointer
          enable(DEFAULT)/disable use of EBP as general purpose register.
          -fno-omit-frame-pointer replaces -fp
-f[no-]exceptions
          enable/disable(DEFAULT) C++ exception handling table generation
-f[no-]exceptions
          enable(DEFAULT)/disable exception handling

Interprocedural Optimization (IPO)
----------------------------------

-[no-]ip  enable(DEFAULT)/disable single-file IP optimization
          within files
-ipo[n]   enable multi-file IP optimization between files
-ipo-c    generate a multi-file object file (ipo_out.o)
-ipo-S    generate a multi-file assembly file (ipo_out.S)
-ip-no-inlining
          disable full and partial inlining
-ip-no-pinlining
          disable partial inlining
-ipo-separate
          create one object file for every source file (overrides -ipo[n])
-ipo-jobs<n>
          specify the number of jobs to be executed simultaneously during the
          IPO link phase

Advanced Optimizations
----------------------

-unroll[n]
          set maximum number of times to unroll loops.  Omit n to use default
          heuristics.  Use n=0 to disable the loop unroller
-[no-]unroll-aggressive
          enables more aggressive unrolling heuristics
-funroll-loops
          unroll loops based on default heuristics
-[no-]scalar-rep
          enable(DEFAULT)/disable scalar replacement (requires -O3)
-[no]pad  enable/disable(DEFAULT) changing variable and array memory layout
-safe-cray-ptr
          Cray pointers do not alias with other variables
-[no-]ansi-alias
          enable/disable(DEFAULT) use of ANSI aliasing rules optimizations;
          user asserts that the program adheres to these rules
-[no-]complex-limited-range
          enable/disable(DEFAULT) the use of the basic algebraic expansions of
          some complex arithmetic operations.  This can allow for some
          performance improvement in programs which use a lot of complex
          arithmetic at the loss of some exponent range.
-reentrancy <keyword>
          specify whether the threaded, reentrant run-time support should be
          used
          Keywords:  none (same as -noreentrancy), threaded, async
-noreentrancy
          do not use threaded, reentrant run-time support
-heap-arrays [n]
          temporary arrays of minimum size n (in kilobytes) are allocated in
          heap memory rather than on the stack.  If n is not specified, 
          all temporary arrays are allocated in heap memory.
-no-heap-arrays
          temporary arrays are allocated on the stack (DEFAULT)
-[no-]opt-multi-version-aggressive
          enables more aggressive multi-versioning to check for pointer
          aliasing and scalar replacement
-opt-ra-region-strategy[=<keyword>]
          select the method that the register allocator uses to partition each
          routine into regions
            routine - one region per routine
            block   - one region per block
            trace   - one region per trace
            loop    - one region per loop
            default - compiler selects best option
-[no-]vec
          enables(DEFAULT)/disables vectorization
-[no-]vec-guard-write
          enables cache/bandwidth optimization for stores under conditionals
          within vector loops
-vec-threshold[n]
          sets a threshold for the vectorization of loops based on the
          probability of profitable execution of the vectorized loop in
          parallel
-opt-malloc-options={0|1|2|3|4}
          specify malloc configuration parameters.  Specifying a non-zero <n>
          value will cause alternate configuration parameters to be set for
          how malloc allocates and frees memory
-opt-jump-tables=<arg>
          control the generation of jump tables
            default - let the compiler decide when a jump table, a series of
                      if-then-else constructs or a combination is generated
            large   - generate jump tables up to a certain pre-defined size
                      (64K entries)
            <n>     - generate jump tables up to <n> in size
          use -no-opt-jump-tables to lower switch statements as chains of
          if-then-else constructs
-fno-jump-tables
          do not generate jump tables for switches and if-then-else statements
-opt-block-factor=<n>
          specify blocking factor for loop blocking
-opt-streaming-stores <arg>
          specifies whether streaming stores are generated
            always - enables generation of streaming stores under the
            assumption that the application is memory bound
            auto   - compiler decides when streaming stores are used (DEFAULT)
            never  - disables generation of streaming stores
-mkl[=<arg>]
          link to the Intel(R) Math Kernel Library (Intel(R) MKL) and bring
          in the associated headers
            parallel   - link using the threaded Intel(R) MKL libraries. This
                         is the default when -mkl is specified
            sequential - link using the non-threaded Intel(R) MKL libraries
            cluster    - link using the Intel(R) MKL Cluster libraries plus
                         the sequential Intel(R) MKL libraries
-[no-]opt-subscript-in-range
          assumes no overflows in the intermediate computation of the
          subscripts

Profile Guided Optimization (PGO)
---------------------------------

-prof-dir <dir>
          specify directory for profiling output files (*.dyn and *.dpi)
-prof-src-root <dir>
          specify project root directory for application source files to
          enable relative path resolution during profile feedback on sources
          below that directory
-prof-src-root-cwd
          specify the current directory as the project root directory for
          application source files to enable relative path resolution during
          profile feedback on sources below that directory
-[no-]prof-src-dir
          specify whether directory names of sources should be 
          considered when looking up profile records within the .dpi file
-prof-file <file>
          specify file name for profiling summary file
-[no-]prof-data-order
          enable/disable(DEFAULT) static data ordering with profiling
-[no-]prof-func-order
          enable/disable(DEFAULT) function ordering with profiling
-[no-]prof-func-groups
          enable(DEFAULT with PGO)/disable function grouping
-prof-gen[=keyword]
          instrument program for profiling.
          Optional keyword may be srcpos or globdata
-no-prof-gen
          disable profiling instrumentation
-prof-use[=<arg>]
          enable use of profiling information during optimization
            weighted  - invokes profmerge with -weighted option to scale data
                        based on run durations
            [no]merge - enable(default)/disable the invocation of the profmerge
                        tool
-no-prof-use
          disable use of profiling information during optimization
-opt-prefetch[=n]
          enable levels of prefetch insertion, where 0 disables.
          n may be 0 through 4 inclusive.  Default is 2.
-no-opt-prefetch
          disable(DEFAULT) prefetch insertion.  Equivalent to -opt-prefetch=0
-p        compile and link for function profiling with UNIX gprof tool
          On IA32 and Intel(r)64, -pg is also valid
-f[no-]instrument-functions
          determine whether function entry and exit points are instrumented
-prof-hotness-threshold=<val>
          set the hotness threshold for function grouping and function ordering
          val indicates percentage of functions to be placed in hot region.
          This option requires -prof-use 
           and -prof-func-groups or -prof-func-order
          

Optimization Reports
--------------------

-vec-report[n]
          control amount of vectorizer diagnostic information
            n=0    no diagnostic information
            n=1    indicate vectorized loops (DEFAULT)
            n=2    indicate vectorized/non-vectorized loops
            n=3    indicate vectorized/non-vectorized loops and prohibiting
                   data dependence information
            n=4    indicate non-vectorized loops
            n=5    indicate non-vectorized loops and prohibiting data
                   dependence information
-opt-report [n]
          generate an optimization report to stderr
            0   disable optimization report output
            1   minimum report output
            2   medium output (DEFAULT when enabled)
            3   maximum report output
-opt-report-file=<file>
          specify the filename for the generated report
-opt-report-phase=<phase>
          specify the phase that reports are generated against
-opt-report-routine=<name>
          reports on routines containing the given name
-opt-report-help
          display the optimization phases available for reporting
-tcheck [mode]
          enable analysis of threaded applications (requires Intel(R) Thread
          Checker; cannot be used with compiler alone)
            tci - instruments a program to perform a thread-count-independent
                  analysis
            tcd - instruments a program to perform a thread-count-dependent
                  analysis (DEFAULT when mode is not used)
            api - instruments a program at the api-imports level
-tprofile
          generate instrumentation to analyze multi-threading performance
          (requires Intel(R) Thread Profiler; cannot be used with compiler
          alone)
-tcollect[=<lib>]
          inserts instrumentation probes calling the Intel(R) Trace Collector
          API.  The library -l<lib> is linked in the default being -lVT
          (requires Intel(R) Trace Collector)
-tcollect-filter file
          Enable or disable the instrumentation of specified functions. 
          (requires Intel(R) Trace Collector)

OpenMP* and Parallel Processing
-------------------------------

-openmp   enable the compiler to generate multi-threaded code based on the
          OpenMP* directives
-openmp-profile
          enable analysis of OpenMP application when the Intel(R) Thread
          Profiler is installed
-openmp-stubs
          enables the user to compile OpenMP programs in sequential mode.  The
          OpenMP directives are ignored and a stub OpenMP library is linked
          (sequential)
-openmp-report{0|1|2}
          control the OpenMP parallelizer diagnostic level
-openmp-lib <ver>
          choose which OpenMP library version to link with
            compat - use the GNU compatible OpenMP run-time libraries
                     (DEFAULT)
-openmp-link <library>
          choose whether to link with the static or dynamic OpenMP
          libraries.  Default is dynamic.
-openmp-threadprivate <ver>
          choose which threadprivate implementation to use
            compat - use the GNU compatible thread local storage
            legacy - use the Intel compatible implementation
                     (DEFAULT)
-cluster-openmp
          allows the user to run an OpenMP program on a cluster. See Cluster
          OMP documentation
-cluster-openmp-profile
          link a Cluster OMP program with profiling information. See Cluster
          OMP documentation
-[no-]clomp-sharable-propagation
          reports variables that need to be made sharable by the user with
          Cluster OMP.  See Cluster OMP documentation
-[no-]clomp-sharable-info
          reports variables that the compiler automatically makes sharable for
          Cluster OMP.  See Cluster OMP documentation
-[no-]clomp-sharable-commons
           makes all COMMONs sharable by default for Cluster OMP.  See Cluster 
           OMP documentation 
-[no-]clomp-sharable-modvars
          makes all variables in modules sharable by default for Cluster OMP.
          See Cluster OMP documentation 
-[no-]clomp-sharable-localsaves
          makes all SAVE variables sharable by default for Cluster OMP.  See
          Cluster OMP documentation 
-[no-]clomp-sharable-argexprs
          makes all expressions in function and subroutine call statements
          sharable by default for Cluster OMP.  See Cluster OMP documentation 
-parallel
          enable the auto-parallelizer to generate multi-threaded code for
          loops that can be safely executed in parallel
-par-report{0|1|2|3}
          control the auto-parallelizer diagnostic level
-par-threshold[n]
          set threshold for the auto-parallelization of loops where n is an
          integer from 0 to 100
-[no-]par-runtime-control
          enable compiler to generate runtime control code for effective
          automatic parallelization
-par-schedule-static[=n]
          Specifies a scheduling algorithm for DO loop iteration.
          Divides iterations into contiguous pieces.  Size n if
          specified, equal sized pieces if not.
-par-schedule-static-balanced[=n]
          Divides iterations into even-sized chunks.  Size n if
          specified, equal sized pieces if not.
-par-schedule-static-steal[=n]
          Divides iterations into even-sized chunks, but allows
          threads to steal parts of chunks from neighboring threads
-par-schedule-dynamic[=n]
          Specifies a scheduling algorithm for DO loop iteration.
          Assigns iterations to threads in chunks dynamically.
          Chunk size is n iterations if specified, otherwise 1.
-par-schedule-guided[=n]
          Specifies a scheduling algorithm for DO loop iteration.
          Indicates a minimum number of iterations.  If specified,
          n is the minimum number, otherwise 1.
-par-schedule-guided-analytical[=n]
          Divides iterations by using exponential distribution or 
          dynamic distributions.
-par-schedule-runtime
          Specifies a scheduling algorithm for DO loop iteration.
          Defers the scheduling decision until runtime.
-par-schedule-auto
          Lets the compiler or run-time system determine the 
          scheduling algorithm.
-par-affinity=[<modifier>,...]<type>[,<permute>][,<offset>] 
          tune application performance by setting different thread affinity
-par-num-threads=<n>
          tune application performance by setting different number of threads 

Floating Point
--------------

-fp-model <name>
          enable <name> floating point model variation
            [no-]except - enable/disable floating point semantics
            fast[=1|2]  - enables more aggressive floating point optimizations
            precise     - allows value-safe optimizations
            source      - enables intermediates in source precision
            strict      - enables -fp-model precise -fp-model except, disables
                          contractions and enables pragma stdc fenv_access
-fp-speculation=<mode>
          enable floating point speculations with the following <mode>
          conditions:
            fast   - speculate floating point operations (DEFAULT)
            safe   - speculate only when safe
            strict - same as off
            off    - disables speculation of floating-point operations
-pc32     set internal FPU precision to 24 bit significand
-pc64     set internal FPU precision to 53 bit significand
-pc80     set internal FPU precision to 64 bit significand (DEFAULT)
-mp1      improve floating-point precision (speed impact less than -mp)
-mieee-fp
          same as -mp, can be disabled with -mno-ieee-fp
-[no-]prec-sqrt
          determine if certain square root optimizations are enabled
-[no-]prec-div
          improve precision of FP divides (some speed impact)
-[no-]fast-transcendentals
          generate a faster version of the transcendental functions
-[no-]fp-port
          round fp results at assignments and casts (some speed impact)
-fp-stack-check
          enable fp stack checking after every function/procedure call
-rcd      rounding mode to enable fast float-to-int conversions
-rounding-mode chopped
          set internal FPU rounding control to truncate
-[no-]ftz
          enable/disable flush denormal results to zero
-fpe{0|1|3}
          specifies program-wide behavior on floating point exceptions
-fpe-all={0|1|3}
          specifies floating point exception behavior on all functions
          and subroutines.  Also sets -assume ieee_fpe_flags
-[no]fltconsistency
          specify that improved floating-point consistency should be used
-[no]recursive
          compile all procedures for possible recursive execution

Inlining
--------

-inline-level=<n>
          control inline expansion:
            n=0  disable inlining 
            n=1  inline functions declared with ATTRIBUTES INLINE or
                   FORCEINLINE
            n=2  inline any function, at the compiler's discretion 
-f[no-]inline-functions
          inline any function at the compiler's discretion
-finline-limit=<n>
          set maximum number of statements a function can have and still be
          considered for inlining
-inline-min-size=<n>
          set size limit for inlining small routines
-no-inline-min-size
          no size limit for inlining small routines
-inline-max-size=<n>
          set size limit for inlining large routines
-no-inline-max-size
          no size limit for inlining large routines
-inline-max-total-size=<n>
          maximum increase in size for inline function expansion
-no-inline-max-total-size
          no size limit for inline function expansion
-inline-max-per-routine=<n>
          maximum number of inline instances in any function
-no-inline-max-per-routine
          no maximum number of inline instances in any function
-inline-max-per-compile=<n>
          maximum number of inline instances in the current compilation
-no-inline-max-per-compile
          no maximum number of inline instances in the current compilation
-inline-factor=<n>
          set inlining upper limits by n percentage
-no-inline-factor
          do not set set inlining upper limits
-inline-forceinline
          treat inline routines as forceinline
-inline-calloc
          directs the compiler to inline calloc() calls as malloc()/memset()

Output, Debug, PCH
------------------

-c        compile to object (.o) only, do not link
-S        compile to assembly (.s) only, do not link
-fsource-asm
          produce assembly file with optional source annotations (requires -S)
-f[no-]verbose-asm
          produce assembly file with compiler comments (DEFAULT) (requires -S)
-fcode-asm
          produce assembly file with optional code annotations (requires -S)
-use-msasm
          support Microsoft* style assembly language insertion using MASM style
          syntax
-o <file>
          name output file
-g        produce symbolic debug information in object file (implies -O0 when
          another optimization option is not explicitly set)
-debug [keyword]
          enable debug information and control output of enhanced debug
          information
            keywords:  all, full, minimal, none, [no]inline-debug-info
                       [no]variable-locations, [no]semantic-stepping,
                       extended
                        parallel
-debug-parameters [keyword]
          control output of debug information for PARAMETERS
            keywords: all, used, none (same as -nodebug-parameters)
-nodebug-parameters
          do not output debug information for PARAMETERS
-g0       disable generation of symbolic debug information
-gdwarf-2
          enable generation of debug information using the DWARF2 format
-[no]d-lines
          compile debug statements (indicated by D in column 1)
-DD       compile debug statements, indicated by D in column 1.  This option
          prevents the definition of a macro named D using the command line
          -Dname option (use -Dname=n syntax instead)
-ftrapuv  trap uninitialized variables
-map-opts
          enable option mapping tool
-print-multi-lib
          print information about libraries being used

Preprocessor
------------

-D<name>[=<text>]
          define macro
-nodefines, -noD
          specifies that any -D macros go to the preprocessor only, and not to
          the compiler
-U<name>  remove predefined macro
-allow nofpp-comments
          If a Fortran end-of-line comment is seen within a #define, treat it
          as part of the definition.  Default is allow:fpp-comments
-E        preprocess to stdout
-EP       preprocess to stdout, omitting #line directives
-P        preprocess to file, omitting #line directives
-preprocess-only
          same as -P
-[no]keep  keep/remove preprocessed file generated by preprocessor as input to
           compiler stage.  Not affected by -save-temps.  Default is -nokeep
-fpp[n], -[no]fpp
           run Fortran preprocessor on source files prior to compilation
             n=0     disable running the preprocessor, equivalent to nofpp
             n=1,2,3 run preprocessor

-cpp[n]   same as -fpp[n]
-module path
           specify path where mod files should be placed and first location to
           look for mod files
-I<dir>   add directory to include file search path
-idirafter<dir>
          add directory to the second include file search path (after -I)
-isystem<dir>
          add directory to the start of the system include path
-X, -nostdinc
          remove standard directories from include file search path
-B<prefix>
          find libraries, headers and executables in <prefix>

Component Control
-----------------

-Qoption,<tool>,<opts>
          pass options <opts> to tool specified by <tool>
-Qlocation,<tool>,<dir>
          set <dir> as the location of tool specified by <tool>
-Qinstall <dir>
          set <dir> as root of compiler installation

Language
--------

-[no]altparam
          specify if alternate form of parameter constant declarations
          (without parenthesis) is recognized. Default is to recognize
          Same as -[no]dps
-assume <keyword>
          specify assumptions made by the optimizer and code generator
          keywords: none, [no]byterecl, [no]buffered_io, 
                    [no]bscc (nobscc same as -nbs), 
                    [no]cc_omp, [no]minus0,
                    [no]dummy_aliases (same as -common-args),
                    [no]ieee_fpe_flags,
                    [no]old_boz, [no]old_logical_ldio,
                    [no]old_maxminloc,
                    [no]old_unit_star, [no]old_xor, 
                    [no]protect_constants, [no]protect_parens, 
                    [no]realloc_lhs, [no]2underscore, 
                    [no]underscore (same as -us),
                    [no]std_mod_proc_name, [no]source_include, 
                    [no]split_common,[no]writeable_strings
-ccdefault <keyword>
          specify default carriage control for units 6 and *
          keywords:  default, fortran, list or none
-[no]check <keyword>
          check run-time conditions.  Default is -nocheck
          keywords: all (same as  -C), none (same as -nocheck), 
                    [no]arg_temp_created, [no]bounds (same as -CB), 
                    [no]format, [no]output_conversion, 
                    [no]pointer (same as -CA), 
                    [no]uninit (same as -CU)
-common-args
          assume "by reference" subprogram arguments may alias one
          another.  Same as -assume dummy_aliases
-e03      issue errors for language elements that are not standard in
          Fortran 2003 (same as -stand f03 -warn stderrors options)
-e95      issue errors for language elements that are not standard in
          Fortran 95 (same as -stand f95 -warn stderrors options)
-e90      issue errors for language elements that are not standard in
          Fortran 90 (same as -stand f90 -warn stderrors options)
-[no]extend-source [<keyword>]
          specify rightmost column for fixed form sources
          keywords: 72 (same as -noextend-source and -72),
                    80 (same as -80),
                   132 (same as -132.  Default if you specify 
                        -extend-source without a keyword.)
          
-fixed    specify source files are in fixed format. Same as -FI 
          -nofixed indicates free format
-free     specify source files are in free format. Same as -FR 
          -nofree indicates fixed format
-[no]mixed-str-len-arg
          indicate whether hidden lengths are passed after their
          character argument or after all arguments.
-names <keyword>
          specify how source code identifiers and external names are
          interpreted.
          keywords:  as_is, lowercase (Default. Same as -lowercase), 
                     uppercase (same as -uppercase)
-[no]pad-source
          make compiler acknowledge blanks at the end of a line
-stand [<keyword>]
          specifies level of conformance with ANSI standard to check
          for.  If keyword is not specified, level of conformance is f03
          keywords: f90 (same as -std90), f95(same as -std95), 
                    f03(same as -std95), none (same as -nostand)
-syntax-only
          perform syntax and semantic checking only (no object file produced)

Compiler Diagnostics
--------------------

-cm     suppress all comment messages, same as -warn nousage
-w        disable all warnings
-W<n>     disable warnings (n = 0) or show warnings (n = 1 DEFAULT, same as
          -warn general)
-w90, -w95
          suppress messages about use of non-standard Fortran
-warn <keyword>
          specifies the level of warning messages issued
            keywords: all, none (same as -nowarn)
                      [no]alignments, [no]declarations, [no]errors,
                      [no]general, [no]ignore_loc, [no]interfaces,
                      [no]stderrors, [no]truncated_source, [no]uncalled,
                      [no]unused, [no]usage
-nowarn   suppress all warning messages
-WB       turn a compile-time bounds check into a warning
-Wcheck   enable more strict diagnostics
-Winline  enable inline diagnostics
-[no]traceback
          specify whether the compiler generates PC correlation data used to
          display a symbolic traceback rather than a hexadecimal traceback at
          runtime failure
-[no]gen-interfaces[:[no]source]
          generate interface blocks for all routines in the file.  Can be
          checked using -warn interfaces
          nosource indicates temporary source files should not be saved
-error-limit <size>
          specify the maximum number of error-level or fatal-level compiler
          errors allowed
-noerror-limit
          set no maximum number on error-level or fatal-level error messages
-diag-enable <v1>[,<v2>,...]
          enable the specified diagnostics or diagnostic groups
-diag-disable <v1>[,<v2>,...]
          disable the specified diagnostics or diagnostic groups
          where <vN> may be individual diagnostic numbers or group names.
          where group names include:
              sc[n]      - perform source code analysis: n=1 for critical
                           errors, n=2 for all errors and n=3 for all errors
                           and warnings
              sc-include - perform source code analysis on include files
              sc-parallel[n] - perform analysis of parallelization in 
                               source code: n=1 for critical errors,
                               n=2 for errors, n=3 for all errors and 
                               warnings
-diag-error <v1>[,<v2>,...]
          output the specified diagnostics or diagnostic groups as errors
-diag-warning <v1>[,<v2>,...]
          output the specified diagnostics or diagnostic groups as warnings
-diag-remark <v1>[,<v2>,...]
          output the the specified diagnostics or diagnostic groups as remarks
-diag-dump
          display the currently enabled diagnostic messages to stdout or to a
          specified diagnostic output file.
-diag-file[=<file>]
          <file> where diagnostics are emitted to.  Not specifying this causes
          messages to be output to stderr
-diag-file-append[=<file>]
          <file> where diagnostics are emitted to. When <file> already exists,
          output is appended to the file
-[no-]diag-id-numbers
          enable(DEFAULT)/disable the diagnostic specifiers to be output in
          numeric form
-diag-error-limit <num>
          specify the maximum number of errors emitted

Miscellaneous
-------------

-[no]logo
          display compiler version information.  /nologo disables the output
-V        display compiler version information
-dumpmachine
          display the target machine only
--version
          display GCC style version information
-[no-]sox
          enable/disable(DEFAULT) saving of compiler options and version in
          the executable
-save-temps
          store the intermediate files in current directory and name them
          based on the source file.  Only saves files that are generated by
          default
-dryrun   show driver tool commands but do not execute tools
-v        show driver tool commands and execute tools
-what     display detailed compiler version information
-watch <keyword>
          tells the driver to output processing information
            keywords: all, none (same as -nowatch), [no]source,
                      [no]cmd (same as -v)
-nowatch  suppress processing information output (DEFAULT)
-Tf<file>
          compile file as Fortran source
-multiple-processes[=<n>]
          create multiple processes that can be used to compile large numbers
          of source files at the same time

Data
----

-i{2|4|8}
          set default KIND of integer and logical variables to 2, 4, or 8
-integer-size <size>
          specifies the default size of integer and logical variables
            size:  16, 32, 64
-r{8|16}  set default size of real to 8 or 16 bytes
-real-size <size>
          specify the size of REAL and COMPLEX declarations, constants,
          functions, and intrinsics
            size: 32, 64, 128
-autodouble
          same as -real-size 64 or -r8
-double-size <size>
          defines the size of DOUBLE PRECISION and DOUBLE COMPLEX declarations,
          constants, functions, and intrinsics
            size:  64, 128
-[no]fpconstant
          extends the precision of single precision constants assigned to
          double precision variables to double precision
-[no]intconstant
          use Fortran 77 semantics, rather than Fortran 90/95, to determine
          kind of integer constants
-auto     make all local variables AUTOMATIC. Same as -automatic
-auto-scalar
          make scalar local variables AUTOMATIC (DEFAULT)
-save     save all variables (static allocation) (same as -noautomatic,
          opposite of -auto)
-[no]zero
          enable/disable(DEFAULT) implicit initialization to zero of local
          scalar variables of intrinsic type INTEGER, REAL, COMPLEX, or
          LOGICAL that are saved and not initialized
-dyncom<common1,common2,...>
          make given common blocks dynamically-allocated
-Zp[n]    specify alignment constraint for structures (n=1,2,4,8,16
          -Zp16 DEFAULT)
-[no]align
          analyze and reorder memory layout for variables and arrays
-align <keyword>
          specify how data items are aligned
            keywords: all (same as -align), none (same as -align),
                      [no]commons, [no]dcommons, [no]records,
                      rec1byte, rec2byte, rec4byte, rec8byte, rec16byte,
                      [no]sequence
-fminshared
          Compilation is for the main executable. Absolute addressing can be
          used and non-position independent code generated for symbols that
          are at least protected
-fstack-security-check
          enable overflow security checks.  
          -f[no-]stack-security-check disables (DEFAULT)
-fstack-protector
          enable stack overflow security checks.
          -f[no-]stack-protector disables (DEFAULT)
-fpic, -fPIC
          generate position independent code (-fno-pic/-fno-PIC is DEFAULT)
-fpie, -fPIE
          generate position independent code that will be linked into an
          executable (-fno-pie/-fno-PIE is DEFAULT)
-[no-]global-hoist
          enable(DEFAULT)/disable external globals are load safe
-f[no-]keep-static-consts
          enable/disable(DEFAULT) emission of static const variables even
          when not referenced
-fpack-struct
          pack structure members together
-f[no-]math-errno
          set ERRNO after calling standard math library functions
-no-bss-init
          disable placement of zero-initialized variables in BSS (use DATA)
-mcmodel=<size>
          use a specific memory model to generate code and store data
            small  - Restricts code and  data to the first 2GB of address 
                     space (DEFAULT)
            medium - Restricts code to the first 2GB; it places no memory
                     restriction on data
            large  - Places no memory restriction on code or data
-convert <keyword>
          specify the format of unformatted files containing numeric data
            keywords: big_endian, cray, ibm, little_endian, native, vaxd, vaxg
-falign-functions=[2|16]
          align the start of functions on a 2 (DEFAULT) or 16 byte boundary
-falign-functions
          align the start of functions to an optimal machine-dependent value.
          -fno-align-functions (DEFAULT) aligns on a 2-byte boundary
-fvisibility=[extern|default|protected|hidden|internal]
          Global symbols (data and functions) will get the visibility 
          attribute given by default. Symbol visibility attributes explicitly
          set in the source code or using the symbol visibility attribute
          file options will override the -fvisibility setting
-fvisibility-extern=<file>
          Space separated symbols listed in the <file> argument will geti
          visibility set to extern
-fvisibility-default=<file>
          Space separated symbols listed in the <file> argument will get
          visibility set to default
-fvisibility-protected=<file>
          Space separated symbols listed in the <file> argument will get
          visibility set to protected
-fvisibility-hidden=<file>
          Space separated symbols listed in the <file> argument will get
          visibility set to hidden
-fvisibility-internal=<file>
          Space separated symbols listed in the <file> argument will get
          visibility set to internal
-fvisibility-inlines-hidden
          mark inline member functions as hidden

Compatibility
-------------

-fpscomp <keyword>
          specify the level of compatibility to adhere to with Fortran
          PowerStation
            keywords: all, none (same as -nofpscomp), [no]filesfromcmd,
                      [no]general, [no]ioformat, [no]ldio_spacing,
                      [no]libs, [no]logicals
-nofpscomp
          no specific level of compatibility with Fortran PowerStation
-1, -onetrip
          execute any DO loop at least once
-f66, -66
          allow extensions that enhance FORTRAN-66 compatibility
-f77rtl   specify that the Fortran 77 specific run-time support should be used
          -nof77rtl disables
-vms      enable VMS I/O statement extensions
-gcc-name=<name>
          name and location of gcc if not where expected
-gxx-name=<name>
          name and location of g++ if not where expected
-gcc-version=<version>
          specify the <version> of gcc compatibility.  Default value matches
          gcc version installed.  Major/Minor versions listed but patch
          levels (i.e. 345) are permissible
            320 - gcc 3.2.x compatibility
            330 - gcc 3.3.x compatibility
            340 - gcc 3.4.x compatibility
            400 - gcc 4.0.x compatibility
            410 - gcc 4.1.x compatibility
            420 - gcc 4.2.x compatibility

Linking/Linker
--------------

-L<dir>   instruct linker to search <dir> for libraries
-l<string>
          instruct the linker to link in the -l<string> library
-shared-intel
          link Intel provided libraries dynamically
-static-intel
          link Intel provided libraries statically
-shared-libgcc
          link libgcc dynamically
-static-libgcc
          link libgcc statically
-dynamic-linker<file>
          select dynamic linker other than the default
-no-cxxlib
          do not link in C++ runtime libraries
-cxxlib[=dir]
          link using C++ run-time libraries provided with gcc dir is an
          optional top-level location for the gcc binaries and libraries
-nodefaultlibs
          do not use standard libraries when linking
-nostartfiles
          do not use standard startup files when linking
-nostdlib
          do not use standard libraries and startup files when linking
-nofor-main
          do not link against Fortran main object.  Used when linking Fortran
          objects with C main program
-static   prevents linking with shared libraries
-shared   produce a shared object
-Bstatic  specify following libraries are linked statically
-Bdynamic
          specify following libraries are linked dynamically
-cxxlib-<mode>
          tell the compiler which C++ run-time libraries to use
            nostd - do not link in standard C++ library
-T <file>
          direct linker to read link commands from <file>
-Xlinker <val>
          pass <val> directly to the linker for processing
-Wa,<o1>[,<o2>,...]
          pass options o1, o2, etc. to the assembler
-Wl,<o1>[,<o2>,...]
          pass options o1, o2, etc. to the linker for processing
-Wp,<o1>[,<o2>,...]
          pass options o1, o2, etc. to the preprocessor
-threads  specify that multi-threaded libraries should be linked against
          -nothreads disables multi-threaded libraries

Linker Specific Options
-----------------------

These options are specific to the linker.  Details can be found in the linker
documentation and man page
-L<dir>
-T<arg>
-h<arg>
-u<arg>
-z<arg>
-i
-r
-s
-N
-Bsymbolic
-Bdynamic
-Bstatic

Deprecated Options
------------------

-Ob                      use -inline-level=<n>
-cxxlib-gcc[=dir]        use -cxxlib[=dir]
-i-dynamic               use -shared-intel
-i-static                use -static-intel
-inline-debug-info       use -debug inline-debug-info
-mp                      use -fp-model <arg>
-use-asm                 no replacement
-A-                      use -U<arg>
-prefetch                use -opt-prefetch
-prof-genx               use -prof-gen=srcpos
-fwritable-strings       use -assume writeable-strings 
-openmp-lib=legacy       use -openmp-lib=compat
-xB                      use -xSSE2 (i32 only)
-xK                      use -xSSE (i32 only)
-axK                     No replacement
-xW                      use -msse2
-axW                     use -msse2
-xN                      use -xSSE2
-axN                     use -axSSE2
-xP                      use -xSSE3
-axP                     use -axSSE3
-xT                      use -xSSSE3
-axT                     use -axSSSE3
-xS                      use -xSSE4.1
-axS                     use -axSSE4.1
-xO                      use -msse3
-func-groups             use -prof-func-groups
-diag-enable sv<n>       use -diag-enable sc<n>
-diag-enable sv-include  use -diag-enable sc-include
-diag-disable sv         use -diag-disable sc
-diag-sv                 use -diag-enable sc
-diag-sv-error           use -diag-disable warning
-diag-sv-include         use -diag-enable sc-include
-diag-sv-level           No replacement
-diag-sv-sup             use -diag-disable <v1>[,<v2>,...]

-help [category]   print full or category help message

Valid categories include
       advanced        - Advanced Optimizations
       codegen         - Code Generation
       compatibility   - Compatibility
       component       - Component Control
       data            - Data
       deprecated      - Deprecated Options
       diagnostics     - Compiler Diagnostics
       float           - Floating Point
       help            - Help
       inline          - Inlining
       ipo             - Interprocedural Optimization (IPO)
       language        - Language
       link            - Linking/Linker
       misc            - Miscellaneous
       opt             - Optimization
       output          - Output
       pgo             - Profile Guided Optimization (PGO)
       preproc         - Preprocessor
       reports         - Optimization Reports
       openmp          - OpenMP and Parallel Processing

Copyright (C) 1985-2009, Intel Corporation.  All rights reserved.
* Other brands and names are the property of their respective owners.

