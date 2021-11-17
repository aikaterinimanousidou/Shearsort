#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void printMatrix(double* M, int n, int row, int col);
void sort(double* a, int N, int i, int d, int type);
void readInput(double* M, int N, FILE *stream_in);
void writeOutput(char* filename, double* M, int N);
int cmpfunc(const void *a, const void *b);

int main(int argc, char *argv[]) {
    if(argc != 3){
        printf("Expected input: shearshort input_file output_file\n");
        exit(0);
    }

    int size, rank, N, chunk, d, *row_index;
    double start_time, stop_time, max_time, *matrix, *local_matrix, *temp_1, *temp_2;

    MPI_Init(&argc, &argv);               // Initialize MPI 
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get my number

    if(rank == 0){
        char* input_filename = argv[1];
    
        // Open inputfile for reading
        FILE *stream_in;
        stream_in = fopen(input_filename, "r");
        if(stream_in == NULL){
            printf("Error: Unable to open file: %s\n", input_filename);
            fclose(stream_in);
            exit(0);
        }

        // Read size of matrix
        fscanf(stream_in, "%d ", &N);

        // Check assumption
        if(N % size != 0){
            printf("ERROR! Invalid size! The matrix size must be divisible with the nr of processes\n");
            exit(0);
        }

        // Read input to matrice
        matrix = (double *)malloc(N*N*sizeof(double));
        readInput(matrix, N, stream_in);

        /*
        // Print initial matrix
        printf("Initial matrix\n");
        printMatrix(matrix, N, N, N);
        */
    }

    // Calculate number of needed iterations
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    d = ceil(log2(N));

    // Calculate local_matrix size and allocate memory
    chunk = N/size;
    local_matrix = (double *)malloc(chunk*N*sizeof(double));
    row_index = (int *)malloc(chunk*sizeof(int));

    // Start timer
	start_time = MPI_Wtime();

    // Determine each ranks original row indeces
    for(int i = 0; i < chunk; i++){
        row_index[i] = rank*chunk + i;
    }

    // Scatter data to local_matrix arrays
    MPI_Scatter(&matrix[0], chunk*N, MPI_DOUBLE, &local_matrix[0], chunk*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Main for-loop
    for(int l = 0; l <= d+1; l++){

        // Sort each row depending on original row index
        for(int i=0; i < chunk; i++){            
            if(row_index[i] % 2 == 0){
                //sort(local_matrix, N, i, d, 0);
                qsort(&local_matrix[i*chunk],chunk,sizeof(double), cmpfunc);
            } else {
                sort(local_matrix, N, i, d, 1);
            }
        }

        // If columns need to be sorted
        if(l <= d){

            // Transpose the local matrix
            temp_1 = (double *)malloc(chunk*N*sizeof(double));
            temp_2 = (double *)malloc(chunk*N*sizeof(double));
            
            // Transpose each local matrix
            for(int i = 0; i < chunk; i++){
                for(int j = 0; j < N; j++){
                    temp_1[j*chunk+i] = local_matrix[i*N+j];
                }
            }
            
            // Re-order the rows
            int block = chunk * chunk;
            for(int i = 0; i < chunk; i++){
                for(int j = 0; j < size; j++){
                    for(int k = 0; k < chunk; k++){
                        temp_2[i*N+j*chunk+k] = temp_1[i*chunk+j*block+k];
                    }
                }
            }         
        
            // Assemble the columns
            for(int i = 0; i < chunk; i++){
                MPI_Alltoall(&temp_2[i*size*chunk], chunk, MPI_DOUBLE, &local_matrix[i*size*chunk], chunk, MPI_DOUBLE, MPI_COMM_WORLD);

            }

            // Sort each column
            for(int i=0; i < chunk; i++){
                //sort(local_matrix, N, i, d, 0);
                qsort(&local_matrix[i*chunk],chunk,sizeof(double), cmpfunc);
            }                  

            // Transpose back each local matrix
            for(int i = 0; i < chunk; i++){
                for(int j = 0; j < N; j++){
                    temp_1[j*chunk+i] = local_matrix[i*N+j];
                }
            }

            // Re-order back the rows
            for(int i = 0; i < chunk; i++){
                for(int j = 0; j < size; j++){
                    for(int k = 0; k < chunk; k++){
                        //temp_2[i*chunk+j*block+k] = temp_1[i*N+j*(size)+k];
                        temp_2[i*N+j*chunk+k] = temp_1[i*chunk+j*block+k];
                    }
                }
            }                       
            
            // Assemble the rows
            for(int i = 0; i < chunk; i++){
                MPI_Alltoall(&temp_2[i*size*chunk], chunk, MPI_DOUBLE, &local_matrix[i*size*chunk], chunk, MPI_DOUBLE, MPI_COMM_WORLD);
            }
        
            // Clean up
            free(temp_1);
            free(temp_2);
        }
    }

    // Gather all local matrices
    MPI_Gather(&local_matrix[0], chunk*N, MPI_DOUBLE, &matrix[0], chunk*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Compute time
	stop_time = MPI_Wtime()-start_time; // stop timer
	MPI_Reduce(&stop_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank == 0){
        
        /*
        // Print final matrix
        printf("Final matrix\n");
        printMatrix(matrix, N, N, N);
        */

        printf("Execution time %f\n", max_time);

        // Write output and clean up
        char* output_filename = argv[2];
        writeOutput(output_filename, matrix, N);
        free(matrix);
    }

    // Clean up
    free(local_matrix);
    free(row_index);

    MPI_Finalize();
}

void sort(double* a, int N, int i, int d, int type){
    if(type == 0){  // normal sort
        for(int k = 0; k < d+1; k++){
            for(int j = 0; j < N-1; j++){
                if(a[i*N+j] > a[i*N+j+1]){
                    double temp=a[i*N+j];
                    a[i*N+j]=a[i*N+j+1];
                    a[i*N+j+1]=temp;
                }
            }
        }
    } else if(type == 1){ // reverse sort  
        for(int k = 0; k < d+1; k++){
            for(int j = 0; j < N-1; j++){
                if(a[i*N+j] < a[i*N+j+1]){
                    double temp=a[i*N+j];
                    a[i*N+j]=a[i*N+j+1];
                    a[i*N+j+1]=temp;
                }
            }
        }
    }
}

void printMatrix(double* M, int n, int row, int col){
  for(int i = 0; i < row; i++){
    for(int j = 0; j < col; j++){
      printf("%10f ", M[i*n+j]);
    }
    printf("\n");
  }
}

void readInput(double* M, int N, FILE *stream_in){ 
  for(int i = 0; i < N*N; i++){
    fscanf(stream_in, "%lf ", &M[i]);
  }
  fclose(stream_in);
}

void writeOutput(char* filename, double* M, int N){
  FILE *stream_out;
  stream_out = fopen(filename,"w");
  if(stream_out == NULL){
    printf("Error: unable to open file: %s\n", filename);
    exit(0);
  }
  for(int i = 0; i < N*N; i++){
    fprintf(stream_out, "%-10f", M[i]);
  }
  fclose(stream_out);
}

 int cmpfunc(const void *a, const void *b)
 {
  if(*(double *)a  <  *(double *)b) return -1;
  if(*(double *)a  == *(double *)b) return 0;
  if(*(double *)a  >  *(double *)b) return 1; 
}