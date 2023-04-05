# TopHap+_v0.1.1  

(Copyright 2023, Authors and Temple University; see license below)
Updated April 4, 2023
==================

## Description
TopHap+ identifies clones and infers clone phylogeny from single-cell sequencing data. See Miura et al. (ref. 1) for the detail. You are free to download, modify, and expand this code under a permissive license similar to the BSD 2-Clause License (see below).

## Dependencies
* Windows
* python 3 (v3.8.10 was tested)
* python packages: 
numpy, 
biopython, 
SciPy, 
pydot

 Note: If the installation of these python packages is not easy, you may want to use Anaconda for Python 3 (https://www.anaconda.com/distribution/). Or you can try python3 -m pip install [package name].
* MEGA (https://www.megasoftware.net/)
* graphviz (https://graphviz.org/)

## How to use TopHapPlus

## 1. Refine single-cell sequences before the application of TopHap+ (-BEAM option)
`python TopHapPlus.py -BEAM input VAF HF`

E.g., 

`python TopHapPlus.py -BEAM [path to BEAMin.meg] 0.01 0.01`

### Input file

The input file is the alignment of observed single-cell sequences with MEGA format. Do not include sequence of a normal cell. Please see Example/BEAMin.meg for an example. 
 
* "T": Mutant allele
* "A": Wild-type allele
* "?": Missing base

If cell sequences do not need to be refined, duplicate the input mega file and name it *_BEAM.meg, e.g., for input BEAMin.meg, the duplicate file name is (e.g., BEAMin_BEAM.meg).  
If cell sequences are refined using another method, name the file as *_BEAM.meg. In this case, a file before the refinement is also necessary. Please use the same format as described above. 
Then,

`python TopHapPlus.py -BEAM [path to BEAMin.meg] 0.01 0.01`

In this way, BEAM refinement step is skiped.

If given cell sequences contain more than one mutation at a given position, please use TopHap (https://github.com/SayakaMiura/TopHap). Please refine the cell sequences before the application of TopHap, if they contain many missing bases and incorrect base assignments. 

## 2. Convert indel matrix to fasta before the application of TopHap+ (-Cas option)
python TopHapPlus.py -Cas input VAF HF VAF2
* VAF2: VAF cut-off to attach mutations through mutation ordering analysis (the number of cells)
E.g., 

`python TopHapPlus.py -Cas [path to Casin.txt] 0.01 0.01 30`

### Input file

The input file is a matrix of indel. Please see Example/Casin.txt for an example. 

* The first column: Cell ID
* The other columns: Sites of indels
* Please assign a number for each indel. 

"0": Without indel (wild-type)

"-": missing 

## output files
Output files are produced in the same directory as the input file. 
* TopHapPlus.nwk

 TopHap+ clone phylogeny
* TopHapPlus.txt

 TopHap+ clone annotation for each cell
* TopHapPlus.png and TopHapPlus.gv

 TopHap+ mutation tree. The number is mutation ID, which is the position in the given input file.

## Reference
[1] Sayaka Miura, Tenzin Dolker, Maxwell Sanderford, and Sudhir Kumar. Improving cellular phylogenies through integrated use of mutation order and optimality principles (2023) Submitted to Cancers

--------
Copyright 2022, Authors and Temple University
BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
