#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spkmeans.h"

int ValidateData(char **, int);
int CalculateNandDim(char *, int *, int *); 
int FileToMatrix(char *, double ***); 

double WeightedEuclideanNorm(double ***, int, int, int);
double MatrixRowSum(double ***, int, int); 
int MatrixMultiplication(double ***, double ***, double ***, int);
int Pivot(double ***, double ***, double ***, int *, int *, int); 
int CheckMatrix(double ***, double ***, int);
int CopyMat(double ***, double ***, int); 
int CompEiganValues(const void*, const void*);
int PrintMatrix(double ***, int, int);

int FindClosestCentroids(double *, double***, int, int);
int checkEpsilon(double ***, double ***, int, int);
double CalcDist(double *, double *, int);

int main(int argc, char *argv[]){

    char* goal;
    int N;
    int Dim;

    double** data_matrix;
    double** weighted_matrix;
    double** diagonal_matrix;
    double** lnorm_matrix;
    double** eigen_vectors_matrix; 
    double* eigan_values;

    int i;

    if (ValidateData(argv, argc) == 0 || CalculateNandDim(argv[2], &N, &Dim) == 0){
        printf("Invalid input");
        return 0;
    }
    data_matrix = AllocateMat(N, Dim);
    
    if (data_matrix == NULL){
        printf("An Error Has Occurred");
        return 0;
    }

    if (FileToMatrix(argv[2], &data_matrix) == 0){
        printf("Invalid Input!");
        return 0;
    }

    goal = argv[1];
    if (strcmp(goal,"wam") == 0){ 
        weighted_matrix = AllocateMat(N, N); 
        WeightedAdjancencyMatrix(&data_matrix, &weighted_matrix, N, Dim);
        PrintMatrix(&weighted_matrix, N, N);
        FreeMat(&data_matrix, N);
        FreeMat(&weighted_matrix, N);
        return 1;
    }
    else if (strcmp(goal,"ddg") == 0){ 
        weighted_matrix = AllocateMat(N, N);
        diagonal_matrix = AllocateMat(N, N);
        WeightedAdjancencyMatrix(&data_matrix, &weighted_matrix, N, Dim);
        DiagonalDegreeMatirx(&weighted_matrix, &diagonal_matrix, N);
        PrintMatrix(&diagonal_matrix, N, N);
        FreeMat(&data_matrix, N);
        FreeMat(&weighted_matrix, N);
        FreeMat(&diagonal_matrix, N);
        return 1;
    }
    else if (strcmp(goal,"lnorm") == 0){ 
        weighted_matrix = AllocateMat(N, N);
        diagonal_matrix = AllocateMat(N, N);
        lnorm_matrix = AllocateMat(N,N);
        WeightedAdjancencyMatrix(&data_matrix, &weighted_matrix, N, Dim);
        DiagonalDegreeMatirx(&weighted_matrix, &diagonal_matrix, N);
        NormalizedGraphLaplasian(&diagonal_matrix, &weighted_matrix, &lnorm_matrix, N);
        PrintMatrix(&lnorm_matrix, N, N);
        FreeMat(&data_matrix, N);
        FreeMat(&weighted_matrix, N);
        FreeMat(&diagonal_matrix, N);
        FreeMat(&lnorm_matrix, N);
        return 1;

    }
    else{ 
        eigen_vectors_matrix = AllocateMat(N, N);
        
        eigan_values = (double *)calloc(N, sizeof(double));
        if (eigan_values == NULL){
            printf("An Error Has Occurred”");
            return 0;
        }
        if (Jacobi(&data_matrix, &eigen_vectors_matrix, &eigan_values, N) == 0){
            printf("An Error Has Occurred”");
            return 0;
        }

        for (i = 0; i <  N-1; i++){
            printf("%.4f,",eigan_values[i]);
        }
        printf("%.4f",eigan_values[N-1]);
        printf("\n");
        PrintMatrix(&eigen_vectors_matrix, N, N);
        FreeMat(&data_matrix, N);
        FreeMat(&eigen_vectors_matrix, N);
        free(eigan_values);
    }
    return 1;
}

/*
* ValidateData: Validate user data input
* -------------------- 
* args:
*   *argv[] : Contains user input as char*
*   argc : Number of elements in argv
* --------------------
*return:
*   1 if valid, else 0
*/
int ValidateData(char *argv[], int argc){
    if (argc != 3 || (strcmp(argv[1],"wam") != 0 && strcmp(argv[1],"ddg") != 0 && strcmp(argv[1],"lnorm") != 0 && strcmp(argv[1],"jacobi") != 0)){
        return 0;
    }
    return 1;
}

/*
* CalculateNandDim: Calculate number of vectors and dimention of vectors
* -------------------- 
* args:
*   char *file_name : The path to the Input file
*   int *N : Pointer to N, the variable will store number of vectors
*   int *Dim : Pointer to Dim, the variable will store the dimention of the vectors
* --------------------
* return:
*   0 if error, else 1
*/
int CalculateNandDim(char *file_name, int *N, int *Dim){
    FILE *ifp = NULL;
    int NumOfVectors = 0;
    int CurrDim = 1;
    char c;

    ifp = fopen(file_name,"r");
    
    if (ifp == NULL){
        return 0;
    }
    while ((c = fgetc(ifp)) != EOF){
        if (c ==',' && NumOfVectors == 0){
            CurrDim++;
        }
        if (c == '\n'){
            NumOfVectors++;
        }
    }

    *N = NumOfVectors;
    *Dim = CurrDim;
    fclose(ifp);
    return 1;
}

/*
* FileToMatrix: Read from Input file and populate matrix
* -------------------- 
* args:
*   char *file_name : The path to the Input file
*   double ***mat : pointer to mat, the variable that will store the matrix
* --------------------
*return:
*   0 if error, else 1
*/
int FileToMatrix(char *file_name, double ***mat){
    FILE *ifp = NULL;
    int i = 0, j = 0, p = 0;
    double cord;
    char c;

    ifp = fopen(file_name,"r");
    if (ifp == NULL){
        return 0;
    }

    while ((EOF != (p = fscanf(ifp, "%lf%c", &cord, &c))) && p == 2) {
        (*mat)[i][j] = cord;
        if (c == ','){
            j++;
        } 
        else{
            j = 0;
            i++;
        }
    }
    fclose(ifp);
    return 1;
}

/*
* AllocateMat: Allocate memory to matrix in size N*Dim
* -------------------- 
* args:
*   int N : Number of rows
*   int Dim : Number of colunms
* --------------------
*return:
*   double **mat : Pointer to the new matix in memory 
*/
double** AllocateMat(int N, int Dim){
    double **mat;
    int i;

    mat = (double **)calloc(N, sizeof(double *));
    if (mat == NULL){
        return NULL;
    }
    for (i = 0; i < N; i++){
        mat[i] = (double*)calloc(Dim, sizeof(double));
        if (mat[i] == NULL){
            return NULL;
        }
    }
    return mat;
}

/*
* FreeMat: free memory of matrix
* -------------------- 
* args:
*   double ***mat : Pointer to matrix
*   int N : Number of rows
* --------------------
*return:
*   0 if error, else 1
*/
int FreeMat(double ***mat, int N){
    int i;
    for (i = 0; i < N; i++){
        free((*mat)[i]);
    }
    free(*mat);
    return 1;
}

/*
* WeightedAdjancencyMatrix: Create and form weighted adjacency matrix from mat
* -------------------- 
* args:
*   double ***mat : Pointer to matrix
*   double ***weighted_matrix : Pointer to new matrix in size of N*N
*   int N : Number of rows in mat
*   int Dim : Number of colunms in mat
* --------------------
*return:
*   1 if succeed
*/
int WeightedAdjancencyMatrix(double ***mat, double ***weighted_matrix, int N, int Dim){
    int i, j;
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if (i == j){
                (*weighted_matrix)[i][j] = 0.0;
            }
            else{
                (*weighted_matrix)[i][j] = WeightedEuclideanNorm(mat,i,j,Dim);
            }
        }
    }
    return 1;
}

/*
* DiagonalDegreeMatirx: Create and form diagonal degree matrix from
*    the weighted adjancency Matrix
* -------------------- 
* args:
*   double ***weighted_matrix : Pointer to the weighted adjancency matrix
*   double ***diagonal_matrix : Pointer to the new matrix in size of N*N 
*   int N : Number of rows in the weighted adjancency Matrix
* --------------------
*return:
*   1 if succeed
*/
int DiagonalDegreeMatirx(double ***weighted_matrix, double ***diagonal_matrix, int N){
    int i, j;
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if (i != j){
                (*diagonal_matrix)[i][j] = 0.0;
            }
            else{
                (*diagonal_matrix)[i][j] = MatrixRowSum(weighted_matrix, i, N);
            }
        }
    }
    return 1;
}

/*
* NormalizedGraphLaplasian: Create The normalized graph Laplacian matrix
* from the weighted adjancency Matrix and the diagonal degree matrix
* -------------------- 
* args:
*   double ***diagonal_matrix : Pointer to the diagonal degree matrix
*   double ***weighted_matrix : Pointer to the weighted adjancency Matrix
*   double ***lnorm_matrix : Pointer to the new matrix in size of N*N,
*   The normalized graph Laplacian matrix before initialization
*   int N : Number of rows in all matrix
* --------------------
*return:
*   1 if succeed
*/
int NormalizedGraphLaplasian(double ***diagonal_matrix, double ***weighted_matrix, double ***lnorm_matrix, int N){
    double ** first_mul;
    int i, j;

    first_mul = AllocateMat(N, N);

    /* Computing D ^ - (0.5) */
    for (i = 0; i < N; i++){
        (*diagonal_matrix)[i][i] = 1 / pow((*diagonal_matrix)[i][i], 0.5);
    }
    
    if (MatrixMultiplication(diagonal_matrix,weighted_matrix,&first_mul, N) == 0){
        return 0;
    }
    
    if (MatrixMultiplication(&first_mul,diagonal_matrix,lnorm_matrix, N) == 0){
        return 0;
    }

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if (i == j){
                (*lnorm_matrix)[i][j] = 1 - ((*lnorm_matrix)[i][j]);
            }
            else{
                (*lnorm_matrix)[i][j] = 0 - ((*lnorm_matrix)[i][j]);
            }
        }
    }
    FreeMat(&first_mul, N);
    return 1;
}

/*
* Eigengap: Implement eigengap heuristic to determine
*   number of clusters k
* -------------------- 
* args:
*   double **eigen_values : List of eigen values
*   int N : Number of eigen values
* --------------------
*return:
*   int index : Number of clusters
*   
*/
int Eigengap(double **eigen_values, int N){
    int i, index = -1;
    double max = -1.0;
    double * eigan_values_copy;
    
    eigan_values_copy = (double *)(calloc(N, sizeof(double)));
    
    if (eigan_values_copy == NULL){
        return -1;
    }

    for (i = 0; i < N; i++){
        eigan_values_copy[i] = (*eigen_values)[i];
    }

    qsort(eigan_values_copy, N, sizeof(double), CompEiganValues);

    for (i = 0; i < floor((N / 2)); i++){
        if (fabs((eigan_values_copy)[i] - (eigan_values_copy)[i+1]) > max){
            max = fabs((eigan_values_copy)[i] - (eigan_values_copy)[i+1]);
            index = i + 1;
        }
    }
    free(eigan_values_copy);
    return index;
}

/*
* CompEiganValues: Compartor for Eigengap function
*   number of clusters k
* -------------------- 
* args:
*   const void *a : first eigen values
*   const void *b : second eigen values
* --------------------
*return:
*   1 if first > second
*   -1 if second > first
*   0 if first == second
*/
int CompEiganValues(const void *a, const void *b){
    double res = *((double *)b) - *((double*)a);
    if (res > 0){
        return 1;
    }
    else{
        if (res < 0){
            return -1;
        }
    }
    return 0;
}

/*
* WeightedEuclideanNorm: Calculate euclidean norm
* -------------------- 
* args:
*   double ***mat : Pointer to matrix
*   int i : Index of the i's row in matrix
*   int j : Index of the j's row in matrix
*   int Dim : Dimentions of vectors
* --------------------
*return:
*   double "": The value in the ij's entery in the weighted adjancency matrix,
*/
double WeightedEuclideanNorm(double ***mat, int i, int j, int Dim){
    double sum = 0;
    int k;

    for (k = 0; k < Dim; k++){
        sum = sum + pow((*mat)[i][k]-(*mat)[j][k],2);
    }
    sum = pow(sum,0.5);
    return exp(-sum / 2);
}

/*
* MatrixRowSum: Sum row in matrix
* -------------------- 
* args:
*   double ***weighted_matrix : Pointer to weighted adjancency matrix
*   int i : Index of row to sum
*   int N : Number of elements in row
* --------------------
*return:
*   double sum: Sum of number in row i 
*/
double MatrixRowSum(double ***weighted_matrix, int i, int N){
    double sum = 0;
    int k;

    for (k = 0; k < N; k++){
        sum = sum + (*weighted_matrix)[i][k];
    }

    return sum;
}

/*
* MatrixMultiplication: Excecute matrix multiplication
* -------------------- 
* args:
*   double ***mat1 : Pointer to first matrix
*   double ***mat2 : Pointer to second matrix
*   double ***res : Pointer to empty matrix that will populate by 
*   the result of marix multiplication 
*   int N : Size of matrixs in N*N
* --------------------
*return:
*   1 if succeed
*/
int MatrixMultiplication(double ***mat1, double ***mat2, double ***res, int N){
    int row1, col1, col2;
    double sum = 0.0;
    
    for (row1 = 0; row1 < N; row1++) {
      for (col2 = 0; col2 < N; col2++) {
        sum = 0.0;
        for (col1 = 0; col1 < N; col1++) {
          sum += (*mat1)[row1][col1] * (*mat2)[col1][col2];
        }
        (*res)[row1][col2] = sum;
      }
    }
    return 1;
}

/*
* Jacobi: Implement jacobi algorithm
* -------------------- 
* args:
*   double ***data_matrix : Symetric real valued square matrix 
*   double ***eigen_vectors_matrix : Output Matrix of the eigan values
*   of the data matrix
*   double **eigan_values : Output List of eigan_values
*   int N :  size os matrixs N*N
* --------------------
*return:
*   1 if succeed, else 0
*/
int Jacobi(double ***data_matrix, double ***eigen_vectors_matrix, double **eigan_values, int N){
    int i, ind1 = 0, ind2 = 0, iter = 0;
    int TOP = 100;
    double ** pivot_matrix, ** temp_matrix, ** A, ** Atag;

    pivot_matrix = AllocateMat(N, N);
    temp_matrix = AllocateMat(N, N);
    A = AllocateMat(N, N);
    Atag = AllocateMat(N, N);

    if (pivot_matrix == NULL || temp_matrix == NULL || A == NULL || Atag == NULL){
        return 0;
    }
    
    CopyMat(&A,data_matrix, N);
    CopyMat(&Atag, &A, N);
    for (i = 0; i < N; i ++){
            pivot_matrix[i][i] = 1.0;
            (*eigen_vectors_matrix)[i][i] = 1.0;
            temp_matrix[i][i] = 1.0;
    }
    
    do {
        if (Pivot(&pivot_matrix, &A, &Atag, &ind1, &ind2, N) == 1){
            MatrixMultiplication(eigen_vectors_matrix, &pivot_matrix, &temp_matrix, N);
            CopyMat(eigen_vectors_matrix, &temp_matrix, N);
            break;
        }
        CopyMat(&A, &Atag, N);
        MatrixMultiplication (eigen_vectors_matrix, &pivot_matrix, &temp_matrix, N);
        CopyMat(eigen_vectors_matrix, &temp_matrix, N);
        iter ++;

    } while (iter < TOP);

    for (i = 0; i < N; i++){
        (*eigan_values)[i] = Atag[i][i];
    }
    FreeMat(&pivot_matrix, N);
    FreeMat(&temp_matrix, N);
    FreeMat(&A, N);
    FreeMat(&Atag, N);
    return 1;
}

/*
* Pivot: Implemet pivot step during jacobi algorithm
* -------------------- 
* args:
*   double ***pivot_matrix : Symetric real valued square matrix 
*   double ***A : Matrix that takes part in calculation
*   double ***Atag : Matrix that takes part in calculation
*   int *ind1 : Index for calculation
*   int *ind2 : Index for calculation
*   int N :  size os matrixs N*N
* --------------------
*return:
*   1 if succeed, else 0
*/
int Pivot(double ***pivot_matrix, double ***A, double ***Atag, int *ind1, int *ind2, int N){
    double max = 0.0;
    double teta = 0.0;
    double t,c,s = 0.0;
    double sign = 1.0;
    int r,k;
    
    /* Reseting the rotation matrix */
    (*pivot_matrix)[*ind1][*ind1] = 1.0;
    (*pivot_matrix)[*ind1][*ind2] = 0.0;
    (*pivot_matrix)[*ind2][*ind1] = 0.0;
    (*pivot_matrix)[*ind2][*ind2] = 1.0;

    /* Finding i,j for the pivot (the largest off-diagonal element) */
    for (r = 0; r < N; r++){
        for (k = r + 1; k < N; k++){
            if (fabs((*A)[r][k]) > max){
                (*ind1) = r;
                (*ind2) = k;
                max = fabs((*A)[r][k]);
            }
        }
    }

    if (max == 0.0){
        return 1;
    }

    /* Calculating the necessary values (teta, t, c, s)*/
    teta = ((*A)[*ind2][*ind2] - (*A)[*ind1][*ind1]) / (2 * (*A)[*ind1][*ind2]);
    sign = teta >= 0 ? 1.0 : -1.0;

    t = sign / (fabs(teta) + pow((teta*teta) + 1.0, 0.5));
    c = 1 / (pow((t*t) + 1.0, 0.5));
    s = t * c;

    (*pivot_matrix)[*ind1][*ind1] = c;
    (*pivot_matrix)[*ind2][*ind2] = c;

    /* Updating the rotation matrix */
    if ((*ind1) > (*ind2)){
        (*pivot_matrix)[*ind1][*ind2] = -s;
        (*pivot_matrix)[*ind2][*ind1] = s;
    }
    else{
        (*pivot_matrix)[*ind2][*ind1] = -s;
        (*pivot_matrix)[*ind1][*ind2] = s;

    }
    /* Calculating Atag matrix */
    for (r = 0; r < N; r++){
        if (r != (*ind1) && r != (*ind2)){
            (*Atag)[r][*ind1] = (c * (*A)[r][*ind1]) - (s * (*A)[r][*ind2]);
            (*Atag)[r][*ind2] = (c * (*A)[r][*ind2]) + (s * (*A)[r][*ind1]);
            (*Atag)[*ind1][r] = (c * (*A)[r][*ind1]) - (s * (*A)[r][*ind2]);
            (*Atag)[*ind2][r] = (c * (*A)[r][*ind2]) + (s * (*A)[r][*ind1]);
        }     
    }
    /* Complete Atag matrix calculation */
    (*Atag)[*ind1][*ind1] = ((c * c) * (*A)[*ind1][*ind1]) + ((s * s) * (*A)[*ind2][*ind2]) - (2 * s * c * (*A)[*ind1][*ind2]);
    (*Atag)[*ind2][*ind2] = ((s * s) * (*A)[*ind1][*ind1]) + ((c * c) * (*A)[*ind2][*ind2]) + (2 * s * c * (*A)[*ind1][*ind2]);
    (*Atag)[*ind1][*ind2] = 0;
    (*Atag)[*ind2][*ind1] = 0;
    return CheckMatrix(A, Atag, N);
}

/*
* CheckMatrix: Check matrix convergens
* -------------------- 
* args:
*   double ***mat1 : Matrix that takes part in convergens calculation
*   double ***mat2 : Matrix that takes part in convergens calculation
*   int N :  size os matrixs N*N
* --------------------
*return:
*   1 if succeed, else 0
*/
int CheckMatrix(double ***mat1, double ***mat2, int N){
    int i, j;
    double temp_sum = 0.0;
    double temp_sum_tag = 0.0 ;
    double eps = 0.00001;

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if (i != j){
                temp_sum += (*mat1)[i][j] * (*mat1)[i][j];
                temp_sum_tag += (*mat2)[i][j] * (*mat2)[i][j];
            }
        }
    }

    return (eps >= (temp_sum - temp_sum_tag));
}

/*
* CopyMat: Copy matrix values to second matrix
* -------------------- 
* args:
*   double ***mat1 : Pointer to first matrix (copied matrix)
*   double ***mat2 : Pointer to second matrix (original matrix)
*   int N : Size of matrixs is N*N
* --------------------
*return:
*   1 if succeed
*/
int CopyMat(double ***mat1, double ***mat2, int N){
    int i, j;

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            (*mat1)[i][j] = (*mat2)[i][j];
        }
    }
    return 1;
}

/*
* SortVectors: Sort vectors in eigan vectors matrix 
* and eigan values 
* -------------------- 
* args:
*   double eigan_vectors_matrix : Eigan vectors matrix
*   double ***mat2 : Eigan values array
*   int N :  size os matrixs N*N
* --------------------
*return:
*   1 if succeed, else 0
*/
int SortVectors(double ***eigan_vectors_matrix, double **eigan_values, int N){
    int i, row, flag = 0;
    double temp;

    while (!flag){
        flag = 1;
        for (i = 0; i < N - 1; i++){
            if ((*eigan_values)[i] < (*eigan_values)[i+1]){
                flag = 0;
                temp = (*eigan_values)[i+1];
                (*eigan_values)[i+1] = (*eigan_values)[i];
                (*eigan_values)[i] = temp;
                for (row = 0; row < N; row++){
                    temp = (*eigan_vectors_matrix)[row][i+1];
                    (*eigan_vectors_matrix)[row][i+1] = (*eigan_vectors_matrix)[row][i];
                    (*eigan_vectors_matrix)[row][i] = temp;
                }
            }
        }
    }
    return 1;
}

/*
* UMatrix: eigan vectors and create U matrix
* from N*N size to N*K
* -------------------- 
* args:
*   double*** eigan_vactors_matrix : Eigan vectors matrix
*   double ***u_matrix : Updated matrix
*   int N :  size of matrixs N*N
*   int K : Updated number of columns
* --------------------
*return:
*   1 if succeed, else 0
*/
int UMatrix(double ***eigan_vactors_matrix, double ***u_matrix, int N, int K){
    int i,j;
    for (i = 0; i < N; i++){
        for (j = 0; j < K; j++){
            (*u_matrix)[i][j] = (*eigan_vactors_matrix)[i][j];
        }
    }
    return 1;
}

/*
* TMatrix: Form the matrix T in size of N*K
*   from U by renormalizing each of U’s rows 
* -------------------- 
* args:
*   double*** u_matrix : Eagan vectors as colunms in size of N*K
*   double *** t_matrix : Updated matrix
*   int N :  size of matrixs N*N
*   int K : Updated number of columns
* --------------------
*return:
*   1 if succeed, else 0
*/
int TMatrix(double *** u_matrix, double *** t_matrix, int N, int K){
    int i, j;
    double row = 0.0;
    
    for (i = 0 ; i < N; i++){
        row = 0.0;
        for (j = 0; j < K; j++){
            row += (*u_matrix)[i][j] * (*u_matrix)[i][j];
        }
        
        row = pow(row, 0.5);
        if (row != 0.0){
            for (j = 0; j < K; j++){
                (*t_matrix)[i][j] = (*u_matrix)[i][j] / row;
            }
        }
    }
    return 1;
}

/*
* PrintMatrix: Print matrix values
* -------------------- 
* args:
*   int rows : Number of rows in matrix
*   int cols : Number of colunms in matrix
* --------------------
*return:
*   1 if succeed
*/
int PrintMatrix(double ***mat, int rows, int cols){
    int row, col;

    for (row = 0; row < rows; row++){
        for (col = 0; col < cols - 1; col++){
            printf("%.4f", (*mat)[row][col]);
            printf(",");
        }
        printf("%.4f", (*mat)[row][cols - 1]);
        printf("\n");
    }
    return 1;
}

/*
* Kmeans: Excecute Kmeans algorithm
* -------------------- 
* args:
*   double ***data_matrix : Pointer to matrix that contains vectors
*   double ***centroids : Pointer to matrix that contains initialized centroids
*   int N : Number of vectors
*   int *Dim : Pointer to Dim, the variable will store the dimention of the vectors
*   int K : Number of clusters
*   int MAX_ITER :
* --------------------
* return:
*   0 if error, else 1
*/
int Kmeans(double ***data_matrix, double ***centroids, int N, int dim, int K, int MAX_ITER){
    int i, j, m, r, iter_counter = 0, index = 0, eps_check = 1;
    double **old_cent, **new_cent;
    int *clusters_size;

    new_cent = AllocateMat(K, dim);
    old_cent = AllocateMat(K, dim);
    clusters_size = (int *)calloc(K, sizeof(int));
    if (new_cent == NULL || old_cent == NULL || clusters_size == NULL){
        return 0;
    }
    /* Initialized first K centroids */
    for (i = 0; i < K; i++){
        for (j = 0; j < dim; j++){
            old_cent[i][j] = (*data_matrix)[i][j];
        }
    }

    while (iter_counter < MAX_ITER && eps_check == 1){
        for (j = 0; j < K; j++){
            for (m = 0; m < dim; m++){
                new_cent[j][m] = 0.0;
            }
        }

        for (j = 0; j < K; j++){
            clusters_size[j] = 0;
        }

        for (i = 0; i < N; i++){
            index = FindClosestCentroids((*data_matrix)[i], &old_cent, K, dim);
            clusters_size[index]++;

            for (r = 0; r < dim; r++){
                new_cent[index][r] += (*data_matrix)[i][r];
            }
        }

        for (i = 0; i < K; i++){
            for(r = 0; r < dim; r++){
                new_cent[i][r] = (new_cent[i][r]) / (clusters_size[i]);
            }
        }
        iter_counter++;
        eps_check = checkEpsilon(&new_cent, &old_cent, K, dim);
        for (i = 0; i < K; i++){
            for (j = 0; j < dim; j++){
                old_cent[i][j] = new_cent[i][j];
            }
        }
    }

    for (i = 0; i < K; i++){
        for (j = 0; j < dim; j++){
            (*centroids)[i][j] = new_cent[i][j];
        }
    }

    FreeMat(&new_cent, K);
    FreeMat(&old_cent, K);
    free(clusters_size);

    return 1;
}

/*
* FindClosestCentroids: Find closest centroid from vector
* -------------------- 
* args:
*   double *vector : Pointer to vector
*   double ***old_cent : Pointer to matrix that contains centroids
*   int K : Number of centroids
*   int dim : Pointer to Dim, the variable will store the dimention of the vectors
* --------------------
* return:
*   int index: Index of closest centroid to *vector in ***old_cent
*/
int FindClosestCentroids(double *vector, double ***old_cent, int k, int dim){
    double temp_dis = -1.0;
    double low_dis = CalcDist(vector, (*old_cent)[0], dim);
    int i, index = 0;

    for (i = 1; i < k; i++){
        temp_dis = CalcDist(vector, (*old_cent)[i], dim);
        if(temp_dis < low_dis){
            index = i;
            low_dis = temp_dis;
        }
    }
    return index;
}

/*
* checkEpsilon: Check if converge condition is centroids achieved
* -------------------- 
* args:
*   double ***old_cent : Pointer to centroids before kmeans iteration
*   double ***new_cent : Pointer to centroids after kmeans iteration
*   int K : Number of centroids
*   int dim : dimention of centroids
* --------------------
* return:
*   1 if converge condition is achieved, else 0
*/
int checkEpsilon(double ***old_cent ,double ***new_cent, int K, int dim){
    double eps = 0.0;
    int i;

    for (i = 0; i < K; i++){
        if (CalcDist((*old_cent)[i], (*new_cent)[i], dim) >= eps){
            return 0;
        }
    }
    return 1;
}

/*
* CalcDist: Calculate euclidian distance between two vectors
* -------------------- 
* args:
*   double *vec1 : first vector
*   double *vec2 : second vector
*   int dim : Dimention of vectors
* --------------------
*return:
*   1 if succeed
*/
double CalcDist(double *vec1, double *vec2, int dim){
    int i;
    double sum = 0.0, temp = 0.0;
    for (i = 0; i < dim; i++){
        temp = vec1[i] - vec2[i];
        temp = temp * temp;
        sum += temp;
    }
    sum = sqrt(sum);
    return sum;
}

