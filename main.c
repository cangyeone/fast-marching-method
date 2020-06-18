#include <unistd.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "heap.h"
#include "fmm.h"

void testheap(void){
    heap h;
    mincoord mr; 
    float key;
    int cond; 
    heap_create(&h, 0, NULL); 
    heap_push(&h, 0.1, mr);
    heap_push(&h, 0.2, mr);
    heap_push(&h, 0.6, mr);
    heap_pop(&h, &key, &mr);
    heap_pop(&h, &key, &mr);
    cond = heap_pop(&h, &key, &mr);
    cond = heap_pop(&h, &key, &mr);
    printf("%d, %f, %d\n", h.active_entries, key, cond);
}

void testhead(void){
    FILE *fp; 
    coord cr;
    float *data;
    fp = fopen("dt", "rb");
    fread(&cr, 4, 10, fp);
    data = (float *)malloc(16*sizeof(float));
    fread(data, 4, 64, fp);
    printf("HEAD:%d,%d,%d,%d\n%f,%f,%f,%f\n", cr.type, cr.n1, cr.n2, cr.n3, cr.b1, cr.b2, cr.b3, cr.d1);
    printf("Data:%f,%f\n%f,%f\n", data[0], data[1], data[2], data[3]);
    fclose(fp);
}
int main(void){
    float b[3];
    fast_marching_main();
}
