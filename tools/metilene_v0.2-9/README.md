# metilene

this is a fork of

[http://legacy.bioinf.uni-leipzig.de/Software/metilene/](http://legacy.bioinf.uni-leipzig.de/Software/metilene/)

## v0.2.9 fixes

- metilene output truncated and unpredictable as of GCC version 12+ due to illegal stack copy operation
- scripts shebang containing absolute path
- pthreads used as linker option only and not as compiler option, old C standard 99
- too small buffer causing potential overflow when pretty printing date and time
- unsafe macro expansions
- tab/space mix of indentations (cosmetic)

### minimal example to reproduce

ensure metilene binary is in `PATH`

```
cd test
./test.sh
```

- v0.2.8 compiled with GCC<=v11 reports 44 DMRs
- v0.2.8 compiled with GCC>=v12 reports 3x DMRs
- v0.2.9 compiled with GCC>=v0 reports 44 DMRs

## if used, please cite

- Jühling F, Kretzmer H, Bernhart SH, Otto C, Stadler PF, Hoffmann S: "metilene: Fast and sensitive calling of differentially methylated regions from bisulfite sequencing data", Genome Res 26.2 (2016) 256-262

## installation

use pre-compiled `metilene.x` or run

```
make
```

## run metilene

- `metilene`: runs the main program
- `metilene_input.pl`: generates an input file for metilene from bed files
- `metilene_output.pl`: postprocesses the output of metilene

See the manual for further instructions.
