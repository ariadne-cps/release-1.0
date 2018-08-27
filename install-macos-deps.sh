#/bin/sh
brew cask uninstall oclint
brew install boost gtk cairo zlib
wget https://sourceforge.net/projects/buddy/files/latest/download
tar xvzf download
cd buddy-2.4
./configure
make
sudo make install
