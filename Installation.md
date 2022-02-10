# INSTALLATION

To "install" the bindigo net program, run

```shell
  . ./install.sh
```

Note that both '.' (dots) are required.
The `install.sh` script  will:
 - build BindoligoNet in the current directory, and
 - updates your `~/.bash_profile` file by appending a line of code to
   that defines and exports the `BINPATH` environment variable. The
   `BINPATH` variable specifies the path where the `bindoligo` programs
   look to read the thermodynamic parameters (i.e., the `bin/` directory
   inside this repository).

If the `bin/` directory is later moved, the `BINPATH` environment
variable will need to be updated.  This is most easily accomplished by
navigating to the new location of the `bin/` directory and running

```shell
  . ./binset.sh
```

(Again, both dots are required.)  Alternatively, if you run:

```shell
  export BINPATH="/full/path/to/bin/"
```
The `BINPATH` variable will be changed within the local shell's environment.

However, since this method of updating the BINPATH variable will only modify the
local shell's environment, these changes will not be permanent.
To make the changes permanent, the `~/.bash_profile` should be updated.
