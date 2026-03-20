#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
int main (){
    const uint64_t N=1<<10;
    float* A=malloc(sizeof(float)*N);
    float* B=malloc(sizeof(float)*N);
    float dot_product= 0.0;
for ( uint64_t i =0 ; i<N; i++){
	A[i]=1.0;
	B[i]=2.0;
}	
for (uint64_t i= 0; i<N; i++){
	dot_product +=A[i] * B[i];
}

printf("El producto escalar es: %.2f/n",dot_product);
free(A);
free(B);
}
