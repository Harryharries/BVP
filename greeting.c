/* Assume a block-data distribution. The rod is of size 1. The rod is partitioned into n sections, each of size k. Time is partitioned into m sections, each of size h. The computations end at time t = 5. The temperature at a given point and time is calculated using the formulas. then print every step for simulation.
 
 Using an SPMD style (not master-slave) of programming,
 develop a parallel version of the problem and implement in MPI.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include <math.h>

#define PI 3.14159265
#define H 0.5
#define K 0.01

double **read_matrix(int row, int col);
void init_matrix(int row, int col, double k, double **mat);
void print_m(int row, int col, double **mat);
void free_matrix(int row, double **mat);

//create the 2d array
double **read_matrix(int row, int col){
    double **mat = (double **) malloc(sizeof(double *)*row);
    int i=0,j=0;
    for(i=0; i<row; i++)
        mat[i] = (double *) malloc(sizeof(double)*col);
    
    for(i=0; i<row; i++){
        for(j=0; j<col; j++){
            mat[i][j] = -1;
        }i
    }
    return mat;
}

//calculate first line for the rod
double *read_array(int col, double k){
    double *array =  (double *) malloc(sizeof(double)*col);
    for (int j = 0; j < col; j++) {
        array[j] = 100 * sin(j * k * PI);
    }
    return array;
}

//print single array
void print_array( int col, double *mat){
    for(int j=0; j<col; j++){
        printf("%lf ",mat[j]); /* Print */
    }
}

//fill the 2d array
void init_matrix(int row, int col,double k, double **mat){
    for(int i = 0; i < 1; i++) {
        for (int j = 0; j < col; j++) {
            mat[i][j] = 100 * sin(j * k * PI);
        }
    }
    for(int i = 1; i<row; i++){
        mat[i][0] = 0;
        mat[i][col-1] = 0;
    }
}

//print the current
void print_m(int row, int col, double **mat){
    for(int i=0; i<row; i++){
        for(int j=0; j<col; j++){
            printf("%lf ",mat[i][j]); /* Print */
        }
        printf("\n");
    }
}

void free_matrix(int row, double **mat){
    int i=0;
    for(i=0;i<row;i++)
        free(mat[i]);
    free(mat);
}

int main(int argc, char* argv[]) {
    int rank, size;     // for storing this process' rank, and the number of processes
    int *sendcounts;    // array describing how many elements to send to each process
    int *displs;        // array describing the displacements where each segment begins
    int sum = 0;                // Sum of counts. Used to calculate displacements
    double rec_buf[1000];          // buffer where the received data should be stored
    
    double start_time; // use these for timing
    double stop_time;
    
    double rod_size = 1;
    double h = H;
    double k = K;
    double t = 5;
    
    double t_size = t/h;
    double r_size = rod_size/k;
    int col = r_size+1; //section_rod_length
    int row = t_size+1; //section_time_length
    
    double r = row/pow(col,2);
    
    double *array;
    array = read_array(col,k);
    
    MPI_Status  status;
    double rec,rec1,rec2,rec3;
    int tag = 0;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int rem = (col)%size;
    sendcounts = malloc(sizeof(int)*size);
    displs = malloc(sizeof(int)*size);
    
    for (int i = 0; i < size; i++) {
        sendcounts[i] = (col)/size;
        if (rem > 0) {
            sendcounts[i]++;
            rem--;
        }
        displs[i] = sum;
        sum += sendcounts[i];
    }
    
    MPI_Scatterv(array, sendcounts, displs, MPI_DOUBLE, &rec_buf,100, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    int mysize = sendcounts[rank];
    double **mat = read_matrix(row,mysize);
    
    
    for (int i = 0; i < mysize; i++){
        mat[0][i] = rec_buf[i];
    }
    
    if(rank == 0){
        // communication start
        start_time = MPI_Wtime();
        for(int i = 0; i < row; i++){
            mat[i][0] = 0;
        }
        
        start_time = MPI_Wtime();
        
        for(int i = 0; i < row; i++){
            for(int j = 0; j < mysize; j++){
                
                if(mat[i][j] == -1){
                    
                    if(j == mysize-1){
                        rec2 = mat[i][j];
                        MPI_Send(&rec2, 1, MPI_DOUBLE, rank+1 , tag, MPI_COMM_WORLD);
                        MPI_Recv(&rec, 1, MPI_DOUBLE, rank+1 ,tag, MPI_COMM_WORLD, &status);
                    }
                    if(j != mysize-1){
                        rec = mat[i-1][j+1];
                    }
                    mat[i][j] = r*mat[i-1][j-1] + (1-2*r)*mat[i-1][j] + r*rec;
                }else if((mat[i][j] == 0) & (j == 0) & (j == mysize-1) & (i !=0 )){
                    rec2 = 0;
                    MPI_Send(&rec2, 1, MPI_DOUBLE, rank+1 , tag, MPI_COMM_WORLD);
                    mat[i][j] = 0;
                }
            }
        }
        
        
    }
    else if(rank == size-1){
        
        for(int i = 0; i < row; i++){
            mat[i][mysize-1] = 0;
        }
        for(int i = 0; i < row; i++) {
            for (int j = 0; j < mysize; j++) {
                
                if(mat[i][j] == -1){
                    if(j == 0){
                        MPI_Recv(&rec, 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD,&status);
                        rec2 = mat[i][j];
                        MPI_Send(&rec2, 1, MPI_DOUBLE, rank-1 , tag, MPI_COMM_WORLD);
                        
                    }
                    if(j != 0){
                        rec = mat[i-1][j-1];
                        
                    }
                    mat[i][j] = r*rec + (1-2*r)*mat[i-1][j] + r*mat[i-1][j+1];
                }else if((mat[i][j] == 0) & (j == 0) & (j == mysize-1) & (i !=0) ){
                    rec2 = 0;
                    MPI_Send(&rec2, 1, MPI_DOUBLE, rank-1 , tag, MPI_COMM_WORLD);
                    mat[i][j] = 0;
                }
            }
        }
        
        
    }
    else{
        for(int i = 0; i < row; i++) {
            for (int j = 0; j < mysize; j++) {
                if(mat[i][j] == -1){
                    if((j == 0) & (j != mysize-1)){
                        MPI_Recv(&rec, 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD,&status);
                        rec1 = mat[i][j];
                        MPI_Send(&rec1, 1, MPI_DOUBLE, rank-1 , tag, MPI_COMM_WORLD);
                        mat[i][j] = r*rec+(1-2*r)*mat[i-1][j]+r*mat[i-1][j+1];
                    }
                    if((j == mysize-1) & (j != 0)){
                        rec2 = mat[i][j];
                        MPI_Send(&rec2, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
                        MPI_Recv(&rec3, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD,&status);
                        mat[i][j] = r*mat[i-1][j-1]+(1-2*r)*mat[i-1][j]+r*rec;
                    }
                    if((j!=0) & (j != mysize-1)){
                        mat[i][j] = r*mat[i-1][j-1]+(1-2*r)*mat[i-1][j]+r*mat[i-1][j+1];
                    }
                    if((j==0) & (j == mysize-1)){
                        MPI_Recv(&rec, 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD,&status);
                        rec1 = mat[i][j];
                        MPI_Send(&rec1, 1, MPI_DOUBLE, rank-1 , tag, MPI_COMM_WORLD);
                        
                        MPI_Send(&rec1, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
                        MPI_Recv(&rec3, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD,&status);
                        
                        mat[i][j] = r*rec+(1-2*r)*mat[i-1][j]+r*rec3;
                    }
                }
            }
        }
    }
    if(rank == 0){
        printf("calculate finish! \n");
        stop_time = MPI_Wtime();
        printf("h= %f\n", h);
        printf("k= %f\n", k);
        printf("q= %d\n", size);
        printf("Total time (sec): %f\n", stop_time - start_time);
        
        double *array;
        printf("the result matrix :\n");
        for (int i = 0; i < row; i++) {
            
            for (int j = 0; j < mysize; j++) {
                printf("%lf ", mat[i][j]);
            }
            for (int s = 1; s<size; s++){
                int curr_size = sendcounts[s];
                array = malloc (curr_size * sizeof(double));
                MPI_Recv(array,curr_size,MPI_DOUBLE,s,tag,MPI_COMM_WORLD,&status);
                print_array(curr_size,array);
                free(array);
            }
            printf("\n");
        }
        
    }else{
        for(int s = 0; s<row; s++){
            MPI_Send(mat[s], mysize, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
        }
        
    }
    
    MPI_Finalize();
    free(array);
    free(sendcounts);
    free(displs);
    return(0);
    
    
}


