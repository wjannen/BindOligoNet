#!/bin/bash

export BINPATH="`pwd | sed '/bin\/\{0,1\}$/!s/\/\{0,1\}$/\/bin\//'`"
sed -i '' '/BINPATH/d' ~/.bash_profile
echo "export BINPATH=$BINPATH" >> ~/.bash_profile

touch ~/.zshrc
sed -i '' '/BINPATH/d' ~/.zshrc
echo "export BINPATH=$BINPATH" >> ~/.zshrc
