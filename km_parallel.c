// kmeans paralelizado desenvolvido por lucas joao e wesley mayk - 2016

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define RANDNUM_W 521288629;
#define RANDNUM_Z 362436069;

//-----------------------------------------------------------------------------
//geracao de nros aleatorios

unsigned int randum_w = RANDNUM_W;
unsigned int randum_z = RANDNUM_Z;

void srandnum(int seed) {
    unsigned int w, z;
    w = (seed * 104623) & 0xffffffff;
    randum_w = (w) ? w : RANDNUM_W;
    z = (seed * 48947) & 0xffffffff;
    randum_z = (z) ? z : RANDNUM_Z;
}

unsigned int randnum(void) {
    unsigned int u;
    randum_z = 36969 * (randum_z & 65535) + (randum_z >> 16);
    randum_w = 18000 * (randum_w & 65535) + (randum_w >> 16);
    u = (randum_z << 16) + randum_w;
    return (u);
}

//fim geracao de nros aleatorios
//-----------------------------------------------------------------------------
//declarao das vars globais
typedef float* vector_t;
int npoints;
int dimension;
int ncentroids;
float mindistance;
int seed;
vector_t *data, *centroids;
int *map;
int *dirty;
int too_far;
int has_changed;
int nro_threads;
int pnts_in_thrds;
int pnts_in_cntrds;
//fim declarao das vars globais
//-----------------------------------------------------------------------------
//calcula distancia entre dois pontos
float v_distance(vector_t a, vector_t b) {
    int i;
    float distance = 0;
    for (i = 0; i < dimension; i++)
        distance +=  pow(a[i] - b[i], 2);
    return sqrt(distance);
}
//fim calcula distancia entre dois pontos
//-----------------------------------------------------------------------------
//popula as particoes
static void *populate(void *arg) {
    int i, j;
    float tmp;
    float distance;

    int arg_thread = (int) (intptr_t) arg;
    int baseCalc = arg_thread * pnts_in_thrds;

    for (i = baseCalc ; i < baseCalc + pnts_in_thrds; i++) {
        distance = v_distance(centroids[map[i]], data[i]);
        /* Look for closest cluster. */
        for (j = 0; j < ncentroids; j++) {
            /* Point is in this cluster. */
            if (j == map[i]) continue;
            tmp = v_distance(centroids[j], data[i]);
            if (tmp < distance) {
                map[i] = j;
                distance = tmp;
                dirty[j] = 1;
            }
        }
        /* Cluster is too far away. */
        if (distance > mindistance)
            too_far = 1;
    }
}

static void populateParalelizado(void) {
    int i;

    pthread_t *threads = (pthread_t *) malloc(sizeof(pthread_t)*(nro_threads));
    too_far = 0;

    for (i = 0; i < nro_threads; i++)
        pthread_create(&threads[i], NULL, populate, (void *) (intptr_t) (i));

    for (i = 0; i < nro_threads; i++)
        pthread_join(threads[i], NULL);

    free(threads);
}
//fim popula as particoes
//-----------------------------------------------------------------------------
//calcula centros das particoes
static void *compute_centroids(void *arg) {
    int i, j, k;
    int population;

    int arg_thread = (int) (intptr_t) arg;
    int baseCalc = arg_thread * pnts_in_cntrds;

    for (i = baseCalc; i < baseCalc + pnts_in_cntrds; i++) {
        if (!dirty[i]) continue;
        memset(centroids[i], 0, sizeof(float) * dimension);
        population = 0;
        for (j = 0; j < npoints; j++) {
            if (map[j] != i) continue;
            for (k = 0; k < dimension; k++)
                centroids[i][k] += data[j][k];
            population++;
        }
        if (population > 1) {
            for (k = 0; k < dimension; k++)
                centroids[i][k] *= 1.0/population;
        }
        has_changed = 1;
    }
}

static void computeCentroidsParalelizado(void){
    has_changed = 0;
    pthread_t *threads = (pthread_t *) malloc(sizeof(pthread_t) * nro_threads);
    int i;
    for (i = 0; i < nro_threads; i++)
        pthread_create(&threads[i], NULL, compute_centroids, (void *) (intptr_t) (i));

    for (i = 0; i < nro_threads; i++)
        pthread_join(threads[i], NULL);

    memset(dirty, 0, ncentroids * sizeof(int));
    free(threads);
}
//fim calcula centros das particoes
//-----------------------------------------------------------------------------
//funcao principal kmeans
int *kmeans(void) {
    int i, j, k;
    too_far = 0;
    has_changed = 0;

    if (!(map  = calloc(npoints, sizeof(int))))
        exit (1);
    if (!(dirty = malloc(ncentroids*sizeof(int))))
        exit (1);
    if (!(centroids = malloc(ncentroids*sizeof(vector_t))))
        exit (1);

    for (i = 0; i < ncentroids; i++)
        centroids[i] = malloc(sizeof(float) * dimension);
    for (i = 0; i < npoints; i++)
        map[i] = -1;
    for (i = 0; i < ncentroids; i++) {
        dirty[i] = 1;
        j = randnum() % npoints;
        for (k = 0; k < dimension; k++)
            centroids[i][k] = data[j][k];
        map[j] = i;
    }

    /* Map unmapped data points. */
    for (i = 0; i < npoints; i++)
        if (map[i] < 0)
            map[i] = randnum() % ncentroids;

    do { /* Cluster data. */
        populateParalelizado();
        computeCentroidsParalelizado();
    } while (too_far && has_changed);

    for (i = 0; i < ncentroids; i++)
        free(centroids[i]);
    free(centroids);
    free(dirty);

    return map;
}
//fim funcao principal kmeans
//-----------------------------------------------------------------------------
//main
int main(int argc, char **argv) {
    int i, j;

    if (argc != 7) {
        printf("Usage: npoints dimension ncentroids mindistance seed threads\n");
        exit (1);
    }

    npoints = atoi(argv[1]);
    dimension = atoi(argv[2]);
    ncentroids = atoi(argv[3]);
    mindistance = atoi(argv[4]);
    seed = atoi(argv[5]);
    nro_threads = atoi(argv[6]);

    pnts_in_thrds = npoints / nro_threads;
    pnts_in_cntrds = ncentroids / nro_threads;

    srandnum(seed);

    if (!(data = malloc(npoints*sizeof(vector_t))))
        exit(1);

    for (i = 0; i < npoints; i++) {
        data[i] = malloc(sizeof(float) * dimension);
        for (j = 0; j < dimension; j++)
            data[i][j] = randnum() & 0xffff;
    }

    map = kmeans();

    for (i = 0; i < ncentroids; i++) {
        printf("\nPartition %d:\n", i);
        for (j = 0; j < npoints; j++)
            if(map[j] == i)
                printf("%d ", j);
    }
    printf("\n");

    free(map);
    for (i = 0; i < npoints; i++)
        free(data[i]);
    free(data);

    return (0);
}
//fim main
//-----------------------------------------------------------------------------
