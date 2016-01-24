Source code for the numerical simulations (including license),
including Julia (http://julialang.org/) scripts to reproduce the
figures appearing in the manuscript.

The figures are plotted using the matplotlib Python library. The Julia
source code files have extension .jl and are ASCII text files (UTF-8)
that can be opened with any text editor. See below for a description
of each file.


BP.jl - main Julia module, defining the functions used in simulations
HH.jl - module for constructing and diagonalizing
        the Harper-Hofstadter Hamiltonian
        
nonabelian_fig.jl   -  Fig. 1
selection_fig.jl    - Figs. 2,3,4
torus_edge_fig.jl   - Figs. 5,6
experimental_fig.jl -  Fig. 7



Copyright (c) 2015-2016 and later, Andrei Berceanu. All rights
reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
