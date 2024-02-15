First, SSH into the SCC

Get into queue and reserve a GPU in an interactive session
qrsh -l gpus=1 -l gpu_c=6.0 -pe omp 4

Change directory to get to demo code (This is where I put it)
cd /projectnb/ec526/admin/brower/OpenACC_demo/Integrate


Load in PGI Compiler & gcc
module load pgi
module load gcc

Compile PGI code integrate.cpp
pgc++ -std=c++11 -Minfo=accel -acc -ta=tesla integrate.cpp -o integrate

Run
./integrate

Compare with Serial

Compile
g++ -O2 serial_integrate.cpp -o serial

Run
./serial


Now we can go on to more NVIDIA demos at (where ever you put it)

cd /projectnb/ec526/admin/brower/OpenACC_demo/Poisson2D

with details describe in the manual

openACC_References/OpenACC_Programming_Guide_0_0.pdf

