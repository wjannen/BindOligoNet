CC = gcc
#CFLAGS = -Wall -Werror
LFLAGS = -lm


# Build the bindoligo and bindoligo-D programs,
# then run the binset script that modifies your ~/.bash_profile
# so that this directory's "bin" folder is added to the BINPATH variable,
# which is then exported
all: bindoligo bindoligod
	. ./bin/binset.sh

bindoligo: bindoligo.c
	$(CC) $(CFLAGS) -o bindoligo bindoligo.c $(LFLAGS)

bindoligod: bindoligo.c
	$(CC) $(CFLAGS) -D DNA -o bindoligo-D bindoligo.c $(LFLAGS)
