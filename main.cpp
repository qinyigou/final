//
//  main.cpp
//  matrix
//
//
#include <cstdio>
#include <ctime>
#include <cstdlib>

#define MAXVALUE (10)
#define MAT_SIZE (64)
#define NUM_SIZE (10)

int* generateRandom(int n) {
    
    int size = n * n;
    int* a = (int*) malloc(size * sizeof(int));
    
    for (int i=0; i < size ; i++) {
        
        a[i] = rand() % MAXVALUE;
    }
    
    return a;
}

// generate sparse matrix with sparse ratio 5%
int* generateSparse(int n) {
    
    int size = n * n;
    int* a = (int*) malloc(size * sizeof(int));
    
    for (int i=0; i < size ; i++) {
        
        if(rand() % 20 == 0) {
            
             a[i] = rand() % MAXVALUE;
        }else{
            
            a[i] = 0;
        }
        
       
    }
    
    return a;
}



// standard naive matrix multiplication
void mult(int* a, int* b, int* c, int n) {
    
    for (int i=0; i < n;  i++) {
        
        int in = i * n;
        
        for(int j=0; j < n; j++) {
            
           int* ele = &c[in + j];
            
            *ele = 0;
            
           for(int k=0; k < n; k++) {
            
                
                *ele += a[in + k] * b[k*n + j];
            }
            
        }
    }
}

void sparseMul(int* a, int* b, int* c, int n) {
    
    for(int i=0; i<n*n ; i++) {
        
        c[i] = 0;
    }
    
    for (int i=0; i < n;  i++) {
        
        int in = i * n;
        
        for(int j=0; j < n; j++) {
            
            int aij = a[in + j];
            
            // if a[i][j] == 0 continue
            if(aij == 0) {
                
                continue;
            }
            
            int jn = j * n;
            
            for(int k=0; k < n; k++) {
                
                
                c[in + k] +=  aij * b[jn + k];
            }
            
        }
    }
    
}

// strassen algorithm
void strassen(int* a, int* b, int* c, int n) {
    
    if(n <= MAT_SIZE) {
        
        mult(a, b, c, n);
        
        return;
    }
    
    int half = n / 2;
    
    int hsize = half * half;
    
    int* x = (int*) malloc(hsize * sizeof(int));
    int* y = (int*) malloc(hsize * sizeof(int));
    
    int k = 0;
    int p = 0;
    int f = 0;
    int h = half * n;

    for (int i=0; i < half; i++) {
        
        f += half;
        h += half;
        
        for (int j=0; j < half; j++) {
            
            x[k] = a[p++];
            
            y[k++] = b[f++] - b[h++];
        }
        
        p+=half;
    }
    
    int* p1 = (int*) malloc(hsize * sizeof(int));
    
    strassen(x, y, p1, half);
    
    p = 0;
    h = half * n;
    k = 0;
    for (int i=0; i < half; i++) {
        
        h += half;
        
        for (int j=0; j < half; j++) {
            
            x[k] = a[p] + b[half + p];
            
            p++;
            
            y[k++] = b[h++];
            
            
        }
        
        p+=half;
    }
    
    int* p2 = (int*) malloc(hsize * sizeof(int));
    strassen(x, y, p2, half);
    
    
    p = 0;
    h = half * n;
    k = 0;
    for (int i=0; i < half; i++) {
        
        for (int j=0; j < half; j++) {
            
            x[k] = a[h] + b[half + h];
            
            h++;
            
            y[k++] = b[p++];
        }
        
        h += half;
        
        p += half;
        
    }
    
    int* p3 = (int*) malloc(hsize * sizeof(int));
    strassen(x, y, p3, half);
    
    
    
    int q = half * n;
    p = 0;
    h = half * n;
    k = 0;
    for (int i=0; i < half; i++) {
        
        h += half;
        
        for (int j=0; j < half; j++) {
            
            x[k] = a[h++];
            
            y[k++] = - b[p++] + b[q++];
        }
        
        p += half;
        q += half;
        
    }
    
    int* p4 = (int*) malloc(hsize * sizeof(int));
    strassen(x, y, p4, half);
    
    
    p = 0;
    h = half * n;
    k = 0;
    for (int i=0; i < half; i++) {
        
        h += half;
        
        for (int j=0; j < half; j++) {
            
            x[k] = a[p] + a[h];
            y[k++] = b[p++] + b[h++];
            
        }
        
        p += half;
        
    }
    
    int* p5 = (int*) malloc(hsize * sizeof(int));
    strassen(x, y, p5, half);
    
    
    p = 0;
    h = half * n;
    k = 0;
    for (int i=0; i < half; i++) {
        
        p += half;
        h += half;
        
        for (int j=0; j < half; j++) {
            
            x[k] = a[p] - a[h];
            y[k++] = b[h - half] + b[h];
            p++;
            h++;
        }
        
    }
    
    int* p6 = (int*) malloc(hsize * sizeof(int));
    strassen(x, y, p6, half);
    
    
    
    p = 0;
    h = half * n;
    k = 0;
    for (int i=0; i < half; i++) {
        
        for (int j=0; j < half; j++) {
            
            x[k] = a[p] - a[h];
            y[k++] = b[p] + b[p + half];
            
            p++;
            h++;
        }
        
        p += half;
        h += half;
        
    }
    
    int* p7 = (int*) malloc(hsize * sizeof(int));
    strassen(x, y, p7, half);
    
    
    k = 0;
    
    int offset = half * n;
    
    for (int i=0; i < half; i++) {
        
        int t = i * n;
        
        for (int j=0; j < half; j++) {
            
            c[t + j] = p5[k] + p4[k] - p2[k] + p6[k];
            
            c[t + j + half] = p1[k] + p2[k];
            
            c[t + j + offset] = p3[k] + p4[k];
            
            c[t + j + offset + half] = p1[k] + p5[k] - p3[k] - p7[k];
            
            k++;
        }
    }
    
    free(p1);
    free(p2);
    free(p3);
    free(p4);
    free(p5);
    free(p6);
    free(p7);
}

// print matrix
void printMat(int* a, int n) {
    
    int k = 0;
    
    for (int i=0; i < n; i++) {
        
        for (int j=0; j < n; j++) {
            
            printf("%d\t", a[k++]);
        }
        
        printf("\n");
    }
    
     printf("\n");
}

void testSparse(int n) {
    
    int* a = generateSparse(n);
    
    int* b = generateRandom(n);
    
    int* c = generateRandom(n);
    
    clock_t start, end;
    
    double interval1, interval2;
    
    start = clock();
    
    mult(a, b, c, n);
    
    end = clock();
    
    interval1 = (end - start) * 1.0 / CLOCKS_PER_SEC;
    
    start = clock();
    
    sparseMul(a, b, c, n);
    
    end = clock();
    
    interval2 = (end - start) * 1.0 / CLOCKS_PER_SEC;
    
    //    printf("n = %d, Run time: %.2f seconds\n", n, interval);
    
    printf( "n = %4d, standard = %.2f, sparse = %.2f\n", n, interval1, interval2);
    
    
    free(a);
    
    free(b);
    
    free(c);
}



// test for matrix size nxn
void testN(int n) {
    
    int* a = generateRandom(n);
    
    int* b = generateRandom(n);
    
    int* c = generateRandom(n);
    
    clock_t start, end;
    
    double interval1, interval2;
    
    start = clock();
    
    mult(a, b, c, n);
    
    end = clock();
    
    interval1 = (end - start) * 1.0 / CLOCKS_PER_SEC;
    
    start = clock();
    
    strassen(a, b, c, n);
    
    end = clock();
    
    interval2 = (end - start) * 1.0 / CLOCKS_PER_SEC;
    
//    printf("n = %d, Run time: %.2f seconds\n", n, interval);
    
    printf( "n = %4d, standard = %.2f, strassen = %.2f\n", n, interval1, interval2);
    
    
    free(a);
    
    free(b);
    
    free(c);
}

// test for various size
void testAllSize()
{
    int n = 2;
    
    for (int i=0; i < NUM_SIZE; i++) {
        
        
        testN(n);
        
        n = n * 2;
    }
    
}

// test for various size for sparse matrix
void testAllSizeSparse()
{
    int n = 2;
    
    for (int i=0; i < NUM_SIZE; i++) {
        
        testSparse(n);
        
        n = n * 2;
    }
    
}


int main(int argc, const char * argv[]) {

    printf("Compare standard with sparse matrix multiplication on sparse matrix\n");
    testAllSizeSparse();
    
    printf("\nCompare standard with strassen matrix multiplication\n");
    testAllSize();
    
    return 0;
}

