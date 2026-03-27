//Diferencias finitas

//Condiciones de contorno

// Condiciones iniciales

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <immintrin.h>

#define L=1.0 //Longitud del dominio

#define T_max=2.0 //Este es el tiempo que simularemos
#define NX=100 // Número de puntos espaciales
#define NT=1000 // Número de pasos en el tiempo
#define c=1.0 //velocidad de la onda
#define PI 3.14159265358979323846

int main(){
// La cuerda que estoy imaginando tiene longitud (L), hay que divir en nodos la cuerda y dejar espacio entre ellos, esto lo llamamos dx es la discretización
double dx=L/(NX-1);
double dt=T_max/NT;
//La información de la onda real no puede viajar a través de la malla más rápido de lo que la computadora la calcula
double c_courant = C_WAVE * (dt / dx); // El número de Courant

double c_sq = c_courant * c_courant;

if (c_courant > 1.0) {
        printf("Error: Inestabilidad numérica. Courant (%.2f) > 1.\n", c_courant);
        return 1;
    }



//Asignación de memoria

// Un puntero es solo una variable que guarda una dirección de memoria (una "etiqueta"), no los datos en sí.
// mm_malloc asegura que el bloque de memoria que se entrega empieza en una dirección que es múltiplo exacto de 32
double *u_prev = (double*)_mm_malloc(NX * sizeof(double), 32);
double *u_curr = (double*)_mm_malloc(NX * sizeof(double), 32);
double *u_next = (double*)_mm_malloc(NX * sizeof(double), 32);

if (u_prev == NULL || u_curr == NULL || u_next == NULL) {
        printf("Error: No hay suficiente memoria RAM.\n");
        return 1;
    }

for (int i=0, i<NX , i++){
double x = i*dx; // La computadora solo entiende el índice i (0, 1, 2, 3...), pero las funciones matemáticas necesitan coordenadas físicas reales. Aquí se multiplica el índice por dx (la distancia entre nodos). Si i = 0, entonces x=0.0 (el inicio de la cuerda).Si i está en la mitad de NX, x valdrá 0.5 (el centro físico de la cuerda).Si i = NX - 1, x valdrá 1.0 (el final de la cuerda, asumiendo que L=1.0).
u_curr[i]=sin(PI*x);
u_prev[i]=u_curr[i]; //Asumimos que la velocidad inicial es 0


}
//Asegurando condiciones de frontera de Dirichlet, (Estas condiciones aplican en los extremos de la onda se les fija un valor)

// al hacer esto aseguramos que en cualquier paso del tiempo , pase lo que pase, los extremos no se van a mover.
u_curr[0] = 0.0; u_curr[NX-1] = 0.0;
u_next[0] = 0.0; u_next[NX-1] = 0.0;
u_prev[0] = 0.0; u_prev[NX-1] = 0.0;


// Cargar constantes en registros AVX
__m256d vec_two  = _mm256_set1_pd(2.0);
__m256d vec_c_sq = _mm256_set1_pd(c_sq);
}
