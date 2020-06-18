#ifndef FMM_H 
#define FMM_H
#define INF_VAL  (1e8)
#include "heap.h"

typedef struct grid {
    int dim;
    int shape[3];
    int pad;
    //float dx[3]; 
} grid; 
typedef struct source {
    int n1;
    int n2; 
    int n3;
    float value;
} source;
typedef struct fmatrix { 
    grid g; 
    float *data;
} fmatrix; 
typedef struct imatrix { 
    grid g; 
    int *data;
} imatrix; 
typedef struct coord {
    int type; 
    int n1, n2, n3; 
    float b1, b2, b3;
    float d1, d2, d3; 
} coord; 
typedef struct fmmdata { 
    int dim;
    fmatrix velo;
    fmatrix T;
    imatrix status; 
    coord cr; 
    heap trail; 
    void (*velo_fn)(fmatrix, fmatrix, imatrix, coord, heap*);
    int *data;
} fmmdata; 
typedef struct header {
    int n1;
    int n2; 
    int n3; 
    float b1; 
    float b2; 
    float b3;
    float d1; 
    float d2; 
    float d3; 
} header;

fmatrix init_fmatrix(grid, float);
imatrix init_imatrix(grid, int);
fmmdata init_fmmdata(fmatrix, coord);
void print_matrix(fmatrix);
void set_source(fmmdata*, source);
void fast_marching(fmmdata*);
void testfmm(void);
void fast_marching_main(void);
void get_velo_coord(fmatrix *velo, coord *cr);
#endif