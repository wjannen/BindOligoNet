# Installation

For installation instructions, please refer to the file
[`Installation.md`](./Installation.md).

# Usage

After building and installing the BindoligoNet programs, the general form
for running them is:;

```shell
  bindoligo <options> [s sequence] [t sequence]
```

or, for DNA-DNA pairings, with:

```shell
  bindoligo-D <options> [s sequence] [t sequence]
```

Both the `s` and `t` sequences can be provided by specifying the path of a file
or by specifying the sequence directly on the command line.
To specify files for `s` and `t`, use switches

```shell
  -s <s_file> -t <t_file>
```

Each file should have one sequence per line and should consist
exclusively of `A`, `C`, `G`, `U`, and `T` characters (case
insensitive).  In all cases, sequences should be entered 5' -> 3'.
Note that sequences entered explicitly on the command line must be
provided following all flags.

BindoligoNet can be run in either single mode (accepting two sequences) or batch
mode.  The mode is specified with the `-m` switch:

```shell
  -m s (single, default)
  -m b (batch)
```

There are two ways of using batch mode: one `s` to many `t`'s,
or *n* `s`'s to *n* `t`'s.
For one `s` with many `t`'s, `s` must be specified on the command line,
and `t` given as a file. For example:

```shell
 -m b -t targets.txt acgugaugaugac
```

In the latter case, both must be specified as files of equal length:

```shell
  -m b -s oligos.txt -t targets.txt
```

BindoligoNet can report optimal alignments, MFEs,
or the free energy landscape dG(j).
Output is selected with the *mandatory* `-o` switch.

```shell
  -o a (optimal alignment)
  -o d (MFE)
  -o l (energy landscape)
```

To facilitate automation, BindoligoNet can output to a file,
specified with the `-f` switch.
If this option is omitted, BindoligoNet defaults to printing to standard output.

```shell
  -f <file> (if omitted, outputs to stdout)
```

BindoligoNet can implement all secondary structure removal penalty schemes
detailed in the [paper](http://rna.williams.edu/tools/BINDIGO2.html)
(dG_5, dG_2, dG_0).  In addition, it is possible to enter
custom 5 parameter or 2 parameter models.

```shell
 -P 5 (dG_5, default)
 -P 2 (dG_2)
 -P 0 (dG_0)
 -P C (custom 5 parameter)
 -P c (custom 2 parameter)
```

If custom 5-parameter (`-P C`) or custom 2-paarmeter (`-P c`) options are used,
parameters are specified using the `-C` flag:

```shell
 -P C -C g0,gA,gC,gG,gU
 -P c -C g0,gL
```

# Examples

```shell
bindoligo -m s -o a ccccccc ggggggggggg
```

will output the optimal alignment between the sequences
s=`ccccccc` and t=`ggggggggggg`.

```shell
bindoligo -m b -o d -t targets.txt uacuuaccuggc
```

will output a list of MFEs from pairing the given oligo (`uacuuaccuggc`)
to each sequences in the file `targets.txt`.
