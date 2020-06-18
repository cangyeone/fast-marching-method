#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "heap.h"
#include "fmm.h"
#include <stdio.h> 
#include <math.h>


#define SET_MAT_VAL2D(M, I, J, V) M.data[(I+M.g.pad)*(M.g.shape[1]+M.g.pad*2)+J+M.g.pad]=V; 
#define GET_MAT_VAL2D(M, I, J) (M.data[(I+M.g.pad)*(M.g.shape[1]+M.g.pad*2)+J+M.g.pad])
#define SET_MAT_VAL3D(M, I, J, K, V) M.data[(I+M.g.pad)*(M.g.shape[1]+M.g.pad*2)*(M.g.shape[2]+M.g.pad*2)+(J+M.g.pad)*(M.g.shape[1]+M.g.pad*2)+K+M.g.pad]=V;
#define GET_MAT_VAL3D(M, I, J, K) M.data[(I+M.g.pad)*(M.g.shape[1]+M.g.pad*2)*(M.g.shape[2]+M.g.pad*2)+(J+M.g.pad)*(M.g.shape[1]+M.g.pad*2)+K+M.g.pad]

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

#define SQUARE(x) ((x)*(x))

//初始化浮点矩阵
fmatrix init_fmatrix(grid g, float init_val){
    fmatrix m;
    m.g = g;
    int pad = g.pad; 
    int dim = g.dim;
    int len=1; 
    for(int i=0;i<dim;i++){
        len = len * (g.shape[i]+pad*2);
    }
    m.data = (float *)malloc(len*sizeof(float));
    for(int i=0;i<len;i++)m.data[i]=init_val;
    for(int i=0;i<g.shape[0];i++){
        for(int j=0;j<g.shape[1];j++){
            if(dim==2){
                m.data[i*(g.shape[1]+pad*2)+j]=init_val;
            }
            if(dim==3){
                for(int k=0;k<g.shape[2];k++){
                m.data[i*(g.shape[2]+pad*2)*(g.shape[1]+pad*2)+j*(g.shape[1]+pad*2)+k]=init_val;
                }
                
            }
        }
    }
    return m; 
}
//初始化整形矩阵
imatrix init_imatrix(grid g, int init_val){
    imatrix m;
    m.g = g;
    int pad = g.pad; 
    int dim = g.dim;
    int len=1; 
    for(int i=0;i<dim;i++){
        len = len * (g.shape[i]+pad*2);
    }
    m.data = (int *)malloc(len*sizeof(int));
    for(int i=0;i<len;i++)m.data[i]=2;
    for(int i=0;i<g.shape[0];i++){
        for(int j=0;j<g.shape[1];j++){
            m.data[(i+pad)*(g.shape[1]+pad*2)+j+pad]=init_val; 
        }
    }
    return m;
}
//打印矩阵
void print_matrix(fmatrix a){
    grid g = a.g;
    int pad = g.pad; 
    int dim = g.dim;
    float val; 
    printf("DIM:%d, %d, %d\n", dim, g.shape[0], g.shape[1]);
    
    for(int i=0;i<g.shape[0];i++){
        for(int j=0;j<g.shape[1];j++){
            val = a.data[(i+pad)*(g.shape[1]+pad*2)+j+pad]; 
            if(val==INF_VAL){
                val = 0.0;
            }
            printf("%5.3f ", val);
        }
        printf("\n");
    }
        
}
//打印矩阵
void print_imatrix(imatrix a){
    grid g = a.g;
    int pad = g.pad; 
    int dim = g.dim;
    printf("DIM:%d, %d, %d\n", dim, g.shape[0], g.shape[1]);
    
    for(int i=0;i<g.shape[0];i++){
        for(int j=0;j<g.shape[1];j++){
            printf("%d ", a.data[(i+pad)*(g.shape[1]+pad*2)+j+pad]);
        }
        printf("\n");
    }
        
}

void velo_2d_0(fmatrix Time, fmatrix velo, imatrix status, coord cr, heap *h){
    float t1, t2; 
    float v, t, crrt;
    float sqr; 
    float dd1, dd2, d12; 
    float d1, d2, d3; 

    d1 = cr.d1;
    d2 = cr.d2; 
    d3 = cr.d3; 
    dd1 = d1 * d1; 
    dd2 = d2 * d2; 
    d12 = d1 * d2;
    v = 1/GET_MAT_VAL2D(velo, cr.n1, cr.n2);
    crrt = GET_MAT_VAL2D(Time, cr.n1, cr.n2);
    t1 = MIN(GET_MAT_VAL2D(Time, cr.n1-1, cr.n2), GET_MAT_VAL2D(Time, cr.n1+1, cr.n2));
    t2 = MIN(GET_MAT_VAL2D(Time, cr.n1, cr.n2-1), GET_MAT_VAL2D(Time, cr.n1, cr.n2+1));

    sqr = dd1 * v * v + dd2 * v * v - (t1 - t2) * (t1 - t2); 
    t = 1 / (dd1+dd2) * (dd1 * t2 + dd2 * t1 + d12 * sqrt(fabs(sqr)));
    if(t > MAX(t1, t2) && sqr > 0){
        t = t;
    }
    else{
        t = MIN(t1+d1*v, t2+d2*v);
    }
    
    SET_MAT_VAL2D(Time, cr.n1, cr.n2, MIN(t, crrt));
    SET_MAT_VAL2D(status, cr.n1, cr.n2, 1); 
    mincoord mc; 
    mc.n1 = cr.n1; mc.n2 = cr.n2; 

    heap_push(h, GET_MAT_VAL2D(Time, mc.n1, mc.n2), mc);
}
float cal_velo_2d(float t1, float t2, float d1, float d2, float v){
    float dd1, dd2, d12; 
    float t, sqr; 
    dd1 = d1 * d1; 
    dd2 = d2 * d2; 
    d12 = d1 * d2;
    sqr = dd1 * v * v + dd2 * v * v - (t1 - t2) * (t1 - t2); 
    t = 1 / (dd1+dd2) * (dd1 * t2 + dd2 * t1 + d12 * sqrt(fabs(sqr)));
    if(t > MAX(t1, t2) && sqr > 0){
        t = t;
    }
    else{
        t = MIN(t1+d1*v, t2+d2*v);
    }
    return t;
}
void velo_3d_0(fmatrix Time, fmatrix velo, imatrix status, coord cr, heap *h){
    float t1, t2, t3; 
    float v, t, crrt;
    float sqr; 
    float dd1, dd2, dd3, d12, d23, d13; 
    float d1, d2, d3; 
    float d123; 
    float Z, d1232, d14, d24, d34; 

    d1 = cr.d1;
    d2 = cr.d2; 
    d3 = cr.d3; 
    dd1 = d1 * d1; 
    dd2 = d2 * d2; 
    dd3 = d3 * d3; 
    d12 = d1 * d2;
    d23 = d2 * d3; 
    d13 = d1 * d3; 
    d123 = d1 * d2 * d3; 
    d14 = dd1 * dd1; 
    d24 = dd2 * dd2; 
    d34 = dd3 * dd3; 
    v = 1/GET_MAT_VAL3D(velo, cr.n1, cr.n2, cr.n3);
    crrt = GET_MAT_VAL3D(Time, cr.n1, cr.n2, cr.n3);
    t1 = MIN(GET_MAT_VAL3D(Time, cr.n1-1, cr.n2, cr.n3), GET_MAT_VAL3D(Time, cr.n1+1, cr.n2, cr.n3));
    t2 = MIN(GET_MAT_VAL3D(Time, cr.n1, cr.n2-1, cr.n3), GET_MAT_VAL3D(Time, cr.n1, cr.n2+1, cr.n3));
    t3 = MIN(GET_MAT_VAL3D(Time, cr.n1, cr.n2, cr.n3-1), GET_MAT_VAL3D(Time, cr.n1, cr.n2, cr.n3+1));
    sqr = dd1 * v * v + dd2 * v * v - (t1 - t2) * (t1 - t2); 
    Z = d12*d12+d23*d23+d13*d13;
    d1232 = d123 * d123;
    if(d1232*v>(d24*dd1*(2*dd1-dd3)*(t1-t3)*(t1-t3))+d34*dd1*(2*dd1-dd2)*(t1-t2)*(t1-t2)){
        t = 1/Z * (d12 * t3 + d13 * t2 + d23 * t1 + d123*sqrt(Z*v-dd3*SQUARE(t1-t2)-dd2*SQUARE(t1-t3)-dd3*SQUARE(t1-t2)));
    }
    else{
        t = cal_velo_2d(t2, t3, d2, d3, v);
    }
    if(d1232*v>(d14*dd2*(2*dd2-dd1)*(t2-t1)*(t2-t1)+d34*dd2*(2*dd2-dd3)*(t2-t3)*(t2-t3))){
        t = 1/Z * (d12 * t3 + d13 * t2 + d23 * t1 + d123*sqrt(Z*v-dd3*SQUARE(t1-t2)-dd2*SQUARE(t1-t3)-dd3*SQUARE(t1-t2)));
    }
    else{
        t = cal_velo_2d(t2, t3, d2, d3, v);
    }
    if(d1232*v>(d14*dd2*(2*dd2-dd1)*(t2-t1)*(t2-t1)+d24*dd3*(2*dd3-dd1)*(t1-t3)*(t1-t3))){
        t = 1/Z * (d12 * t3 + d13 * t2 + d23 * t1 + d123*sqrt(Z*v-dd3*SQUARE(t1-t2)-dd2*SQUARE(t1-t3)-dd3*SQUARE(t1-t2)));
    }
    else{
        t = cal_velo_2d(t2, t3, d2, d3, v);
    }


    SET_MAT_VAL3D(Time, cr.n1, cr.n2, cr.n3, MIN(t, crrt));
    SET_MAT_VAL3D(status, cr.n1, cr.n2, cr.n3, 1); 
    mincoord mc; 
    mc.n1 = cr.n1; mc.n2 = cr.n2; mc.n3 = cr.n3;

    heap_push(h, GET_MAT_VAL3D(Time, mc.n1, mc.n2, mc.n3), mc);
}
//初始化数据
fmmdata init_fmmdata(fmatrix velo, coord cr){
    grid g = velo.g;
    fmmdata data; 
    heap h; 
    data.velo = velo;
    
    //printf("SSS");
    data.T = init_fmatrix(g, INF_VAL);
    data.status = init_imatrix(g, 0);
    //printf("DDDI:%d, %d\n", data.T.g.shape[0], data.T.g.shape[1]);
    data.dim = velo.g.dim;
    if(data.dim==2){
        if(data.cr.type==0){
            data.velo_fn = velo_2d_0;
        }
    }
    if(data.dim==3){
        if(data.cr.type==0){
            data.velo_fn = velo_3d_0;
        }
    }    
    
    data.cr = cr; 
    
    heap_create(&h, 0, NULL);
    data.trail = h;
    /*
    return data; */
    return data;
}


//设置源
void set_source(fmmdata *fd, source a){
    if(fd->dim==2){
        coord cr;
        float val; 
        SET_MAT_VAL2D(fd->status, a.n1, a.n2, 2); 
        SET_MAT_VAL2D(fd->T, a.n1, a.n2, a.value);
        
        cr = fd->cr; 
        cr.n1 = a.n1;
        cr.n2 = a.n2-1; 
        fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
        cr = fd->cr; 
        cr.n1 = a.n1;
        cr.n2 = a.n2+1; 
        fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
        cr = fd->cr; 
        cr.n1 = a.n1-1;
        cr.n2 = a.n2; 
        fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
        cr = fd->cr; 
        cr.n1 = a.n1+1;
        cr.n2 = a.n2; 
        fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
    }
    else if(fd->dim==3){
        coord cr;
        float val; 
        SET_MAT_VAL3D(fd->status, a.n1, a.n2, a.n3, 2); 
        SET_MAT_VAL3D(fd->T, a.n1, a.n2, a.n3, a.value);
        
        cr = fd->cr; 
        cr.n1 = a.n1;
        cr.n2 = a.n2-1; 
        cr.n3 = a.n3;
        fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
        cr = fd->cr; 
        cr.n1 = a.n1;
        cr.n2 = a.n2+1; 
        cr.n3 = a.n3;
        fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
        cr = fd->cr; 
        cr.n1 = a.n1-1;
        cr.n2 = a.n2; 
        cr.n3 = a.n3;
        fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
        cr = fd->cr; 
        cr.n1 = a.n1+1;
        cr.n2 = a.n2; 
        cr.n3 = a.n3;
        fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
        cr = fd->cr; 
        cr.n1 = a.n1;
        cr.n2 = a.n2; 
        cr.n3 = a.n3+1;
        fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
        cr = fd->cr; 
        cr.n1 = a.n1+1;
        cr.n2 = a.n2; 
        cr.n3 = a.n3-1;
        fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
    }
}

void fast_marching(fmmdata *fd){
    coord cr;
    float key; 
    mincoord mr; 
    float val; 
    int cond; 
    int count = 0;
    if(fd->dim==2){
        while(1){
            cond = heap_pop(&fd->trail, &key, &mr);
            //break;
            if(cond==0)break; 
            if(GET_MAT_VAL2D(fd->status, mr.n1, mr.n2)==2){
                continue;
            }
            else{
                SET_MAT_VAL2D(fd->status, mr.n1, mr.n2, 2);
                cr = fd->cr;
                cr.n1 = mr.n1;
                cr.n2 = mr.n2-1;
                if(GET_MAT_VAL2D(fd->status, cr.n1, cr.n2)!=2)fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
                cr.n1 = mr.n1;
                cr.n2 = mr.n2+1;
                if(GET_MAT_VAL2D(fd->status, cr.n1, cr.n2)!=2)fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
                cr.n1 = mr.n1-1;
                cr.n2 = mr.n2;
                if(GET_MAT_VAL2D(fd->status, cr.n1, cr.n2)!=2)fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
                cr.n1 = mr.n1+1;
                cr.n2 = mr.n2;
                if(GET_MAT_VAL2D(fd->status, cr.n1, cr.n2)!=2)fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
            }

        }
    }
    else if(fd->dim==3){
        while(1){
            cond = heap_pop(&fd->trail, &key, &mr);
            //break;
            if(cond==0)break; 
            if(GET_MAT_VAL3D(fd->status, mr.n1, mr.n2, mr.n3)==2){
                continue;
            }
            else{
                SET_MAT_VAL3D(fd->status, mr.n1, mr.n2, mr.n3, 2);
                cr = fd->cr;
                cr.n1 = mr.n1;
                cr.n2 = mr.n2-1;
                cr.n3 = mr.n3;
                if(GET_MAT_VAL3D(fd->status, cr.n1, cr.n2, cr.n3)!=2)fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
                cr.n1 = mr.n1;
                cr.n2 = mr.n2+1;
                cr.n3 = mr.n3;
                if(GET_MAT_VAL3D(fd->status, cr.n1, cr.n2, cr.n3)!=2)fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
                cr.n1 = mr.n1-1;
                cr.n2 = mr.n2;
                cr.n3 = mr.n3;
                if(GET_MAT_VAL3D(fd->status, cr.n1, cr.n2, cr.n3)!=2)fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
                cr.n1 = mr.n1+1;
                cr.n2 = mr.n2;
                cr.n3 = mr.n3;
                if(GET_MAT_VAL3D(fd->status, cr.n1, cr.n2, cr.n3)!=2)fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
                cr.n1 = mr.n1;
                cr.n2 = mr.n2;
                cr.n3 = mr.n3-1;
                if(GET_MAT_VAL3D(fd->status, cr.n1, cr.n2, cr.n3)!=2)fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
                cr.n1 = mr.n1;
                cr.n2 = mr.n2;
                cr.n3 = mr.n3+1;
                if(GET_MAT_VAL3D(fd->status, cr.n1, cr.n2, cr.n3)!=2)fd->velo_fn(fd->T, fd->velo, fd->status, cr, &fd->trail);
            }

        }
    }
}


void get_velo_coord(fmatrix *velo, coord *cr){
    FILE *fp; 
    float *data;
    int dim=3;
    int len=1;
    fmatrix vlo; 
    grid g; 
    fp = fopen("dt", "rb");
    fread(cr, 4, 10, fp);
    if(cr->n3==0){
        dim=2;
    }
    g.shape[0] = cr->n1;
    g.shape[1] = cr->n2; 
    g.shape[2] = cr->n3;
    g.pad = 2;
    g.dim = dim;
    for(int i=0;i<dim;i++){
        len = len * (g.shape[i]);
    }
    vlo = init_fmatrix(g, 1.0);
    //printf("value\n");
    data = (float *)malloc(len*sizeof(float));
    fread(data, sizeof(float), len, fp);
    for(int i=0;i<g.shape[0];i++){
        for(int j=0;j<g.shape[1];j++){
            if(dim==2){
                SET_MAT_VAL2D(vlo, i, j, data[i*g.shape[1]+j]);
            }
            else{
                for(int k=0;k<g.shape[2];k++){
                    SET_MAT_VAL3D(vlo, i, j, k, data[i*g.shape[2]*g.shape[1]+j*g.shape[1]+k]);
                }
            }
        }
    }
    *velo = vlo;
    velo->g=g; 
    free(data);
    data = NULL; 
    fclose(fp);
}
void get_source(fmmdata *fmm){
    FILE *fp; 
    int num;
    source sr;
    fp = fopen("ds", "rb");
    fread(&num, 4, 1, fp);
    printf("NumS:%d\n", num);
    for(int i=0;i<num;i++){
        fread(&sr, sizeof(source), 1, fp);
        set_source(fmm, sr);
    }
    fclose(fp);
    //print_matrix(fmm.T);
    //print_imatrix(fmm->status);
}
void save_matrix(fmatrix m, coord cr){
    grid g=m.g; 
    int len=1; 
    float *data;
    FILE *fp; 
    int dim=m.g.dim;
    for(int i=0;i<g.dim;i++){
        len = len * (g.shape[i]);
    }
    fp = fopen("d-io", "wb");
    
    //fwrite(&cr, 4, 1, fp);
    data = (float *)malloc(len*sizeof(float));
    fread(data, sizeof(float), len, fp);
    for(int i=0;i<g.shape[0];i++){
        for(int j=0;j<g.shape[1];j++){
            if(g.dim==2){
                data[i*g.shape[1]+j]=GET_MAT_VAL2D(m, i, j);
            }
            else{
                for(int k=0;k<g.shape[2];k++){
                    data[i*g.shape[2]*g.shape[1]+j*g.shape[1]+k]=GET_MAT_VAL3D(m, i, j, k);
                }
            }
        }
    }
    if(dim==2)cr.n3=0;

    fwrite(data, sizeof(float), len, fp);
    free(data);
    data=NULL;
    fclose(fp);
}
void fast_marching_main(void){
    grid g;
    fmatrix velo;
    fmmdata data; 
    coord cr;
    source sr; 
    fmmdata fmm; 
    source s;
    printf("Get velo\n");
    get_velo_coord(&velo, &cr);
    fmm = init_fmmdata(velo, cr);
    get_source(&fmm);
    fast_marching(&fmm);
    //print_matrix(fmm.T);
    save_matrix(fmm.T, cr);
}

void testfmm(void){
    grid g;
    fmatrix velo;
    fmmdata data; 
    coord cr;
    source sr; 
    fmmdata fmm; 
    source s;
    g.dim = 2; 
    g.shape[0] = g.shape[1] = g.shape[2] = 8;
    g.pad = 2;
    velo = init_fmatrix(g, 1.0); 
    cr.d1 = cr.d2 = cr.d3 = 1.; 
    cr.d1 = 0.1; cr.d2 = 0.2;
    cr.n1 = cr.n2 = cr.n3 = 8;
    cr.type=0;
    printf("BBB%d\n", fmm.dim);
    fmm = init_fmmdata(velo, cr);
    printf("AAA%d\n", fmm.dim);
    sr.n1 = 1; 
    sr.n2 = 1; 
    sr.n3 = 1;
    sr.value = 0.; 
    printf("AAA%d\n", fmm.dim);
    set_source(&fmm, sr);
    //print_matrix(data.velo);
    fast_marching(&fmm);
    print_matrix(fmm.T);
    print_imatrix(fmm.status);

}