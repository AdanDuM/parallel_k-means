#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <pthread.h>
#include <stdint.h>

/*!
 *  começo geracao nros aleatorios (todos os pontos)
 */
#define RANDNUM_W 521288629;
#define RANDNUM_Z 362436069;

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
// fim geracao nros aleatorios (todos os pontos)

/*!
 *  comeco declaracao de vars
 */
typedef float* vector_t;

int npoints;
int dimension;
int ncentroids;               //!< nro de particoes
float mindistance;
int seed;                     //!< semente utilizada para gerar nros
// data possui os pontos
// centroids possui o centro da particao
vector_t *data, *centroids;
int *map;                     //!< associa cada ponto a uma particao
int *dirty;                   //!< define se particao esta suja ou limpa
int too_far;
int has_changed;
int nro_threads;
int pnts_in_thrds;
int pnts_in_cntrds;
pthread_mutex_t mutex_compute_centroids;
pthread_mutex_t mutex_parcial_kmeans;
// fim declarao de vars

/*!
 *  comeco calculo kmeans
 */
float v_distance(vector_t a, vector_t b) { //!< calcula dist. entre dois pontos
  int i;
  float distance = 0;
  for (i = 0; i < dimension; i++)
    distance +=  pow(a[i] - b[i], 2);
  return sqrt(distance);
}

static void populate(int arg_thread) {
  int i, j;
  float tmp;
  float distance;
  too_far = 0;
  // associa cada ponto a cada centro de particao
  for (i = arg_thread * pnts_in_thrds; i < (arg_thread * pnts_in_thrds) + pnts_in_thrds; i++) {
    // pthread_mutex_lock(&mutex_parcial_kmeans);
      distance = v_distance(centroids[map[i]], data[i]);
    // pthread_mutex_unlock(&mutex_parcial_kmeans);
    for (j = arg_thread * pnts_in_cntrds; j < (arg_thread * pnts_in_cntrds) + pnts_in_cntrds; j++) {
    // for (j = 0; j < ncentroids; j++) {
      // so executa se o ponto nao for daquela particao
      if (j == map[i]) continue;
      tmp = v_distance(centroids[j], data[i]);
      // se uma distancia melhor, entao muda particao do ponto
      if (tmp < distance) {
        map[i] = j;
        distance = tmp;
        dirty[j] = 1;
      }
    }
    // verifica se clusterizacao aceitavel
    if (distance > mindistance)
      too_far = 1;
  }
}

static void compute_centroids(void) {
  pthread_mutex_lock(&mutex_compute_centroids);
  int i, j, k;
  int population;
  has_changed = 0;
  for (i = 0; i < ncentroids; i++) {
    // so executa se particao estiver suja
    if (!dirty[i]) continue;
    // zera centro das particoes
    memset(centroids[i], 0, sizeof(float) * dimension);
    // calcula centro da particao
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
  memset(dirty, 0, ncentroids * sizeof(int)); //!< todas particoes limpas
  pthread_mutex_unlock(&mutex_compute_centroids);
}

// void startMemKmeans(void) {
//   int i;

//   // aloca memoria
//   if (!(map  = calloc(npoints, sizeof(int))))
//     exit (1);
//   if (!(dirty = malloc(ncentroids*sizeof(int))))
//     exit (1);
//   if (!(centroids = malloc(ncentroids*sizeof(vector_t))))
//     exit (1);

//   for (i = 0; i < ncentroids; i++)
//     centroids[i] = malloc(sizeof(float) * dimension);
// }

void doMapAndCentroids() {
  int i, j, k;

  for (i = 0; i < npoints; i++)
    map[i] = -1;                      //!< todos pontos nao mapeados

  for (i = 0; i < ncentroids; i++) {
    dirty[i] = 1;                     //!< particoes estao "sujas"
    j = randnum() % npoints;
    for (k = 0; k < dimension; k++)
      centroids[i][k] = data[j][k];   //!< def pontos centros de particoes
    map[j] = i;
  }
}

void *kmeans(void *arg) {
  int i, j, k;
  int arg_thread = (int) (intptr_t) arg;
  too_far = 0;
  has_changed = 0;

  // for (i = arg_thread * pnts_in_thrds; i < (arg_thread * pnts_in_thrds) + pnts_in_thrds; i++)
  //   map[i] = -1;                      //!< todos pontos nao mapeados

  // pthread_mutex_lock(&mutex_parcial_kmeans);
  //   for (i = 0; i < ncentroids; i++) {
  //     dirty[i] = 1;                     //!< particoes estao "sujas"
  //     j = randnum() % npoints;
  //     for (k = 0; k < dimension; k++)
  //       centroids[i][k] = data[j][k];   //!< def pontos centros de particoes
  //     map[j] = i;
  //   }
  // pthread_mutex_unlock(&mutex_parcial_kmeans);

  // pontos nao mapeados recebem particao
  for (i = arg_thread * pnts_in_thrds; i < (arg_thread * pnts_in_thrds) + pnts_in_thrds; i++)
    if (map[i] < 0)
      map[i] = randnum() % ncentroids;

  do { // realiza kmeans
    populate(arg_thread);
    // pthread_mutex_lock(&mutex_compute_centroids);
      compute_centroids();
    // pthread_mutex_unlock(&mutex_compute_centroids);
  } while (too_far && has_changed);

}
// fim calculo kmeans

void finishMemKmeans(void) {
  int i;

  // libera memoria
  for (i = 0; i < ncentroids; i++)
    free(centroids[i]);
  free(centroids);
  free(dirty);
}

int main(int argc, char **argv) {
  int i, j, tmp;

  // trabalha com args
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

  pthread_t thread[nro_threads];
  pnts_in_thrds = npoints / nro_threads;
  pnts_in_cntrds = ncentroids / nro_threads;

  // gera matriz de dados com pontos
  srandnum(seed);

  if (!(data = malloc(npoints*sizeof(vector_t))))
    exit(1);

  for (i = 0; i < npoints; i++) {
    data[i] = malloc(sizeof(float) * dimension);
    for (j = 0; j < dimension; j++)
      data[i][j] = randnum() & 0xffff;
  }

  // int i;

  // aloca memoria
  if (!(map  = calloc(npoints, sizeof(int))))
    exit (1);
  if (!(dirty = malloc(ncentroids*sizeof(int))))
    exit (1);
  if (!(centroids = malloc(ncentroids*sizeof(vector_t))))
    exit (1);

  for (i = 0; i < ncentroids; i++)
    centroids[i] = malloc(sizeof(float) * dimension);

  doMapAndCentroids();

  // realiza o kmeans
  pthread_mutex_init(&mutex_parcial_kmeans, NULL);
  pthread_mutex_init(&mutex_compute_centroids, NULL);

  // startMemKmeans();

  for (i = 0; i < nro_threads; i++)
    pthread_create(&thread[i], NULL, kmeans, (void *) (intptr_t) (i));

  for (i = 0; i < nro_threads; i++)
    pthread_join(thread[i], NULL);

  finishMemKmeans();

  pthread_mutex_destroy(&mutex_parcial_kmeans);
  pthread_mutex_destroy(&mutex_compute_centroids);

  // printa resultado na tela
  for (i = 0; i < ncentroids; i++) {
    printf("\nPartition %d:\n", i);
    for (j = 0; j < npoints; j++)
      if(map[j] == i)
        printf("%d ", j);
  }
  printf("\n");

  // limpa memoria
  free(map);
  for (i = 0; i < npoints; i++)
    free(data[i]);
  free(data);

  return (0);
}
