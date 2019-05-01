git clone https://github.com/allisonvacanti/XtalComp
cd XtalComp
mkdir build
cd build
cmake ..
make
cd ../..

xtalcomp_root=$(pwd)/XtalComp
xtalcomp=$xtalcomp_root/build/libXtalComp.a
spglib=$xtalcomp_root/build/spglib/libspglib.a
g++ -I$xtalcomp_root main.cpp $xtalcomp $spglib -o compcell

rm -rf XtalComp
