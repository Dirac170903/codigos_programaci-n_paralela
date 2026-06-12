#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <immintrin.h>


#define L 1.0        // Longitud del dominio
#define T_max 2.0    // Tiempo a simular
#define NX 100       // Número de puntos espaciales
#define NT 1000      // Número de pasos en el tiempo
#define C_WAVE 1.0   // Velocidad de la onda unificada
#define PI 3.14159265358979323846

int main() {
double dx = L / (NX - 1.0);
double dt = T_max / NT;
    
double c_courant = C_WAVE * (dt / dx); // Fórmula de la condición de Courant-Friedrichs-Lewy  
double c_sq = c_courant * c_courant; // Al cuadrado, asi lo calculamos una vez evita hacer la multiplicación repetidamente dentro de los bucles

if (c_courant > 1.0) {
  printf("Inestabilidad numérica. Courant (%.2f) > 1.\n", c_courant);
  return 1;
}

double *u_prev = (double*)_mm_malloc(NX * sizeof(double), 32);
double *u_curr = (double*)_mm_malloc(NX * sizeof(double), 32);
double *u_next = (double*)_mm_malloc(NX * sizeof(double), 32);

if (u_prev == NULL || u_curr == NULL || u_next == NULL) {
   printf("Error: No hay suficiente memoria RAM.\n");
   return 1;
}

for (int i = 0; i < NX; i++) {
  double x = i * dx; 
  u_curr[i] = sin(PI * x);
  u_prev[i] = u_curr[i]; 
}
// La onda comienza con la forma de una curva seno. Se asume que la onda parte del reposo, por lo que el estado pasado es igual al presente.


// Condiciones de frontera
u_curr[0] = 0.0; u_curr[NX-1] = 0.0;
u_next[0] = 0.0; u_next[NX-1] = 0.0;
u_prev[0] = 0.0; u_prev[NX-1] = 0.0;

 // Registros AVX. Se preparan registros AVX de 256bits, cada uno almacena 4 valores double simultáneamente.
__m256d vec_two  = _mm256_set1_pd(2.0);
__m256d vec_c_sq = _mm256_set1_pd(c_sq);

//
// BUCLE DE EVOLUCIÓN TEMPORAL
// 
    
for (int t = 0; t < NT; t++) {
        int i = 1;
        
        // Bucle espacial vectorizado (AVX)
        for (; i <= NX - 5; i += 4) {
            // Cargar datos no alineados (ya que i empieza en 1)
            __m256d uc = _mm256_loadu_pd(&u_curr[i]);
            __m256d up = _mm256_loadu_pd(&u_prev[i]);
            __m256d uc_left = _mm256_loadu_pd(&u_curr[i-1]);
            __m256d uc_right = _mm256_loadu_pd(&u_curr[i+1]);

            // Operaciones matemáticas de la ec. de onda
            __m256d term1 = _mm256_sub_pd(_mm256_mul_pd(vec_two, uc), up);
            __m256d term2 = _mm256_add_pd(uc_left, uc_right);
            term2 = _mm256_sub_pd(term2, _mm256_mul_pd(vec_two, uc));
            term2 = _mm256_mul_pd(vec_c_sq, term2);

            __m256d un = _mm256_add_pd(term1, term2);
            
            // Guardar el resultado
            _mm256_storeu_pd(&u_next[i], un);
        }

        // Bucle escalar para el residuo (si los puntos interiores no son múltiplo de 4)
        for (; i < NX - 1; i++) {
            u_next[i] = 2.0 * u_curr[i] - u_prev[i] + c_sq * (u_curr[i-1] - 2.0 * u_curr[i] + u_curr[i+1]);
        }

        // Intercambio de punteros (para avanzar en el tiempo)
        double *temp = u_prev;
        u_prev = u_curr;
        u_curr = u_next;
        u_next = temp;
    }

    printf("Simulación completada con éxito.\n");

    // Liberar memoria alineada
    _mm_free(u_prev);
    _mm_free(u_curr);
    _mm_free(u_next);

    return 0;
}
