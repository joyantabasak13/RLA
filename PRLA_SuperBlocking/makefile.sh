# /opt/homebrew/Cellar/open-mpi/4.1.5/bin/mpicxx -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /opt/homebrew/Cellar/open-mpi/4.1.5/include/  -o hello ./openMPI_Test.cpp
# /opt/homebrew/Cellar/open-mpi/4.1.5/bin/mpirun -np 2 ./hello

# g++ -Wall -o main main.cpp -static
# /opt/homebrew/Cellar/open-mpi/4.1.5/bin/mpicxx -std=c++17 -Wall -O3 -I /usr/local/boost_1_80_0/ -I /opt/homebrew/Cellar/open-mpi/4.1.5/include/  -o helloStatic ./openMPI_Test.cpp -static -v


# /opt/homebrew/Cellar/open-mpi/4.1.5/bin/mpicxx -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /opt/homebrew/Cellar/open-mpi/4.1.5/include/ -o openMPI_rla ./openMPI_rla.cpp
# /opt/homebrew/Cellar/open-mpi/4.1.5/bin/mpirun -np 3 ./openMPI_rla ds1_50k

#  g++ -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /usr/local/  -o eColiGenomeComp 

/opt/homebrew/Cellar/open-mpi/4.1.5/bin/mpicxx -std=c++17 -O3 -I /usr/local/boost_1_80_0/ -I /opt/homebrew/Cellar/open-mpi/4.1.5/include/ -o distributed_RadixSort ./distributed_RadixSort.cpp
/opt/homebrew/Cellar/open-mpi/4.1.5/bin/mpirun -np 1 ./distributed_RadixSort ds11_5M_fl
