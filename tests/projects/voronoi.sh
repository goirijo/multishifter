g++ --std=c++17 -L/home/mesto/programming/CASMcode-dev/build/04x/lib -I/home/mesto/programming/CASMcode-dev/build/04x/include voronoi.cpp -lcasm
LD_LIBRARY_PATH=/home/mesto/programming/CASMcode-dev/build/04x/lib ./a.out

# python3 voronoi.py
# exit()

# frame=1000
# for i in {0..3}; do
#     echo $i

#     python3 ./voronoi.py $frame &
#     let frame=frame+1

#     python3 ./voronoi.py $frame &
#     let frame=frame+1
# done
