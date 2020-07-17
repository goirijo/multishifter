# g++ --std=c++17 -L/home/mesto/programming/CASMcode-dev/build/04x/lib -I/home/mesto/programming/CASMcode-dev/build/04x/include voronoi.cpp -lcasm
# LD_LIBRARY_PATH=/home/mesto/programming/CASMcode-dev/build/04x/lib ./a.out $1
# python3 voronoi.py

frame=1000
for i in {0..359}; do
    echo $i

    LD_LIBRARY_PATH=/home/mesto/programming/CASMcode-dev/build/04x/lib ./a.out $i.2
    python3 ./voronoi.py
    mv ./tmp.png ./vorimation/img$frame.png
    let frame=frame+1

    LD_LIBRARY_PATH=/home/mesto/programming/CASMcode-dev/build/04x/lib ./a.out $i.7
    python3 ./voronoi.py
    mv ./tmp.png ./vorimation/img$frame.png
    let frame=frame+1
done
