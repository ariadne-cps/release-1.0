#/bin/sh
export PKG_CONFIG_PATH=/usr/local/opt/libffi/lib/pkgconfig 
brew install lcov boost cairo zlib
wget https://sourceforge.net/projects/buddy/files/latest/download
tar xvzf download
cd buddy-2.4
./configure
make
sudo make install
