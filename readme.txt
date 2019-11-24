This program uses MPI for parallel computing for better efficiency.
So before you run it, please installing open_MPI from https://www.open-mpi.org/


Then you can:

compile:
        make
run:
        qsub myjob




If you want to change the value h and k
please go to greeting.c

(or you can modify the code and change it to manually input): 

in the c file,
there is one line code on the top:

#define H 0.25
#define K 0.001

you can replce the number with any number you want.



change process number:
please go to myjob

there is one line code :

#$ -pe mpi 8

you can replece the 8 with any number you want.