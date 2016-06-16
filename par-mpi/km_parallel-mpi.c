// kmeans paralelizado com mpi desenvolvido por lucas joao e wesley mayk - 2016

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

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
int ncentroids;               //!< nro de particoes
float mindistance;
int seed;                     //!< semente utilizada para gerar nros
/* data possui os pontos */
/* centroids possui o centro da particao */
vector_t *data, *centroids;
int *map;                     //!< associa cada ponto a uma particao
int *dirty;                   //!< define se particao esta suja ou limpa
int too_far;
int has_changed;
int size;
int rank;
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
static void populate(void) {
  int i, j;
  float tmp;
  float distance;
  too_far = 0;

  int baseCalc = rank * npoints / size;
  int finalCalc = (rank + 1) * npoints / size;

  /* associa cada ponto a cada centro de particao */
  for (i = baseCalc; i < finalCalc; i++) {
    distance = v_distance(centroids[map[i]], data[i]);
    for (j = 0; j < ncentroids; j++) {
      /* so executa se o ponto nao for daquela particao */
      if (j == map[i]) continue;
      tmp = v_distance(centroids[j], data[i]);
      /* se uma distancia melhor, entao muda particao do ponto */
      if (tmp < distance) {
        map[i] = j;
        distance = tmp;
        dirty[j] = 1;
      }
    }
    /* verifica se clusterizacao aceitavel */
    if (distance > mindistance)
      too_far = 1; //!< too_far PRECISA SER ENVIADO
  }
}
//fim popula as particoes
//-----------------------------------------------------------------------------
//calcula centros das particoes
static void compute_centroids(void) {
  int i, j, k;
  int population;
  has_changed = 0;

  int baseCalc = rank * npoints / size;
  int finalCalc = (rank + 1) * npoints / size;

  for (i = 0; i < ncentroids; i++) {
    /* so executa se particao estiver suja */
    if (!dirty[i]) continue;
    /* zera centro das particoes */
    memset(centroids[i], 0, sizeof(float) * dimension);
    /* calcula centro da particao */
    population = 0;
    for (j = baseCalc; j < finalCalc; j++) {
      if (map[j] != i) continue;
      for (k = 0; k < dimension; k++)
        centroids[i][k] += data[j][k];
      population++;
    }
    if (population > 1) {
      for (k = 0; k < dimension; k++)
        centroids[i][k] *= 1.0/population; // centroids PRECISA SER ENVIADO
    }
    has_changed = 1; //!< has_changed PRECISA SER ENVIADO
  }
  /* todas particoes limpas */
  memset(dirty, 0, ncentroids * sizeof(int));
}
//fim calcula centros das particoes
//-----------------------------------------------------------------------------
//funcao principal kmeans
int* kmeans(void) {
  int i, j, k, l;
  too_far = 0;
  has_changed = 0;

  if (!(map  = calloc(npoints, sizeof(int)))) {
    MPI_Finalize();
    exit (1);
  }
  if (!(dirty = malloc(ncentroids*sizeof(int)))) {
    MPI_Finalize();
    exit (1);
  }
  if (!(centroids = malloc(ncentroids*sizeof(vector_t)))) {
    MPI_Finalize();
    exit (1);
  }

  for (i = 0; i < ncentroids; i++)
    centroids[i] = malloc(sizeof(float) * dimension);
  for (i = 0; i < npoints; i++)
    map[i] = -1;                      //!< todos pontos nao mapeados
  for (i = 0; i < ncentroids; i++) {
    dirty[i] = 1;                     //!< particoes estao "sujas"
    j = randnum() % npoints;
    for (k = 0; k < dimension; k++)
      centroids[i][k] = data[j][k];   //!< def pontos centros de particoes
    map[j] = i;
  }

  for (i = 0; i < npoints; i++)       //!< pontos nao mapeados recebem particao
    if (map[i] < 0)
      map[i] = randnum() % ncentroids;

  do {
    populate();
    MPI_Barrier(MPI_COMM_WORLD);
    compute_centroids();
    MPI_Barrier(MPI_COMM_WORLD);
  } while (too_far && has_changed);

  if (rank != 0) {
    MPI_Send(&map, npoints, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
  } else {
    int *map_tmp;
    int rank_tmp;

    if (!(map_tmp  = calloc(npoints, sizeof(int)))) {
      MPI_Finalize();
      exit (1);
    }

    for (i = 0; i < size-1; i++) {
      MPI_Recv(&map_tmp, npoints, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&rank_tmp, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      int baseCalc = rank_tmp * npoints / size;
      int finalCalc = (rank_tmp + 1) * npoints / size;

      for (l = baseCalc; l < finalCalc; l++)
        map[l] = map_tmp[l];
    }

    free(map_tmp);
  }

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
  int i, j, tmp;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // trabalha com args
  if (argc != 6) {
    if (rank == 0)
      printf("Usage: npoints dimension ncentroids mindistance seed\n");
    MPI_Finalize();
    exit(1);
  }

  npoints = atoi(argv[1]);
  dimension = atoi(argv[2]);
  ncentroids = atoi(argv[3]);
  mindistance = atoi(argv[4]);
  seed = atoi(argv[5]);

  srandnum(seed);

  if (!(data = malloc(npoints*sizeof(vector_t)))) {
    MPI_Finalize();
    exit(1);
  }

  for (i = 0; i < npoints; i++) {
    data[i] = malloc(sizeof(float) * dimension);
    for (j = 0; j < dimension; j++)
      data[i][j] = randnum() & 0xffff;
  }

  map = kmeans();

  if (rank == 0) {
    for (i = 0; i < ncentroids; i++) {
      printf("\nPartition %d:\n", i);
      for (j = 0; j < npoints; j++)
        if(map[j] == i)
          printf("%d ", j);
    }
    printf("\n");
  }

  free(map);
  for (i = 0; i < npoints; i++)
    free(data[i]);
  free(data);

  MPI_Finalize();
  return (0);
}
//fim main
//-----------------------------------------------------------------------------
