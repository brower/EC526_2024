First, SSH into the SCC

Running in GPUs

https://www.bu.edu/tech/support/research/software-and-programming/programming/multiprocessor/gpu-computing/

INTERACTIVE vs BATCH
https://www.bu.edu/tech/support/research/system-usage/running-jobs/

INTERACTIVE: 
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

cd /projectnb/ec526/admin/brower/OpenACC/Poisson2D

must of course first load modules 
Here there are scripts  to
(0) Makefile to compile
(1) to ask for GPU: get_gpu_node.sh
(2) Run the interative code:  submit_gpu.sh

v0 runs with 1.14 speed up
v1 speed up 180.27
v2 speed up 23.00
v3 speed up 176.69


with details describe in the manual

openACC_References/OpenACC_Programming_Guide_0_0.pdf


BATCH STILL TRYING TO FIGURE THIS OUT


pgc++ -std=c++11 -Minfo=accel -acc -ta=tesla integrate.cpp -o integrate


qsub -l gpus=1 -b y mycode

or 

qsub -l gpus=1 -l gpu_c=7.0 -pe omp 8 your_batch_script

CHECKING ON JOB

https://www.bu.edu/tech/support/research/system-usage/running-jobs/tracking-jobs/
-l g
qstat -u userID
