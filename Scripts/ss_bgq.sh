# Performance test
qsub -A MPCS51087 -t 5 -n 32 --mode c64 ./nbody-parallel.exe -n 102400 -i 10 -t 1
qsub -A MPCS51087 -t 5 -n 32 --mode c1 ./nbody-parallel.exe -n 102400 -i 10 -t 64

# Strong scaling
qsub -A MPCS51087 -t 8 -n 32   --mode c1 ./nbody-parallel.exe -n 524288 -i 10 -t 64
qsub -A MPCS51087 -t 5 -n 64   --mode c1 ./nbody-parallel.exe -n 524288 -i 10 -t 64
qsub -A MPCS51087 -t 5 -n 128  --mode c1 ./nbody-parallel.exe -n 524288 -i 10 -t 64
qsub -A MPCS51087 -t 5 -n 256  --mode c1 ./nbody-parallel.exe -n 524288 -i 10 -t 64
qsub -A MPCS51087 -t 5 -n 512  --mode c1 ./nbody-parallel.exe -n 524288 -i 10 -t 64
qsub -A MPCS51087 -t 5 -n 1024 --mode c1 ./nbody-parallel.exe -n 524288 -i 10 -t 64

qsub -A MPCS51087 -t 8 -n 32   --mode c64 ./nbody-parallel.exe -n 524288 -i 10 -t 1
qsub -A MPCS51087 -t 5 -n 64   --mode c64 ./nbody-parallel.exe -n 524288 -i 10 -t 1
qsub -A MPCS51087 -t 5 -n 128  --mode c64 ./nbody-parallel.exe -n 524288 -i 10 -t 1
qsub -A MPCS51087 -t 5 -n 256  --mode c64 ./nbody-parallel.exe -n 524288 -i 10 -t 1
qsub -A MPCS51087 -t 5 -n 512  --mode c64 ./nbody-parallel.exe -n 524288 -i 10 -t 1
qsub -A MPCS51087 -t 5 -n 1024 --mode c64 ./nbody-parallel.exe -n 524288 -i 10 -t 1

# Production test
qsub -A MPCS51087 -t 30 -n 1024 --mode c1 ./nbody-parallel.exe -n 1048576 -i 10 -t 64

# Production Simulation
qsub -A MPCS51087 -t 60 -n 1024 --mode c1 ./nbody-parallel.exe -n 1048576 -i 400 -t 64
