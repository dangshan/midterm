//
//  main.cpp
//  FFT2.0
//
//  Created by Nancy Wood on 2015/5/6.
//  Copyright (c) 2015年 Nancy Wood. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int FFT(double*x_r,double*x_i,double*y_r,double*y_i,int N);

int main()

    // insert code here...
    
    {
        // y_k = sum(x_n * w^{-kn}, n=0..N-1)
        // w = cos(2*pi/N)+isin(2*pi/N))
        
        int k, n, N, p;
        double *y_r, *y_i, *x_r, *x_i, w_r, w_i, t1,t2;
        
        printf("Please input p=");
        scanf("%d", &p);
        N = 1 << p;
        printf("N=%d\n",N);
        
        // "malloc" create memory for x and y
        x_r = (double *) malloc(N*sizeof(double));
        x_i = (double *) malloc(N*sizeof(double));
        y_r = (double *) malloc(N*sizeof(double));
        y_i = (double *) malloc(N*sizeof(double));
        
        /*
        for(n=0;n<N;++n)
        {
            x_r[n] = n;
            x_i[n] = 0;
        }
        t1 = clock();
        for(k=0;k<N;++k)
        {
            y_r[k] = 0.0;
            y_i[k] = 0.0;
            for(n=0;n<N;n++)
            {
                w_r = cos(-k*n*2*M_PI/N);
                w_i = sin(-k*n*2*M_PI/N);
                y_r[k] = y_r[k] + x_r[n]*w_r - x_i[n]* w_i;
                y_i[k] = y_i[k] + x_r[n]*w_i + x_i[n]* w_r;
            }
        }
        t2 = clock();
        printf("%f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
        for(n=0;n<N;++n)
        {
            //printf("%d : %f + %f i\n", n, y_r[n], y_i[n]);     剛開始跑來檢查答案，之後數字大就別印了...
        }
        

        system("pause");
        */
        for(n=0;n<N;++n)
        {
            x_r[n] = n;
            x_i[n] = 0;
        }
        
        t1 = clock();
        FFT(x_r,x_i,y_r,y_i,N);
        t2 = clock();
        printf("%f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
        for(n=0;n<N;++n)
        {
            //printf("%d : %f + %f i\n", n, y_r[n], y_i[n]);
        }
        
        printf("Yeah~~~~~\n");
        return 0;
    }
    
    int FFT(double*x_r,double*x_i,double*y_r,double*y_i,int N)
    {
        if(N==1)
        {
            y_r[0] = x_r[0];
            y_i[0] = x_i[0];
            return 0;
        }
        // x = x_r + i * x_i  原本的FT
        // y = y_r + i * x_i  重排整理過的
        
        int k, n;
        // double *u_r, *u_i, *v_r, *v_i, w_r, w_i;
        double *even_r, *even_i, *odd_r, *odd_i, w_r, w_i;
        double *even_FT_r, *even_FT_i, *odd_FT_r, *odd_FT_i;
        
        even_r = (double *) malloc(N*sizeof(double));
        even_i = (double *) malloc(N*sizeof(double));
        odd_r = even_r + N/2; // (double *) malloc(N/2*sizeof(double));
        odd_i = even_i + N/2; //(double *) malloc(N/2*sizeof(double));
        even_FT_r = (double *) malloc(N*sizeof(double));
        even_FT_i = (double *) malloc(N*sizeof(double));
        odd_FT_r = even_FT_r + N/2; //(double *) malloc(N/2*sizeof(double));
        odd_FT_i = even_FT_i + N/2; //(double *) malloc(N/2*sizeof(double));
        
        for(n=0;n<N/2;n++)
        {
            even_r[n] = x_r[2*n];
            even_i[n] = x_i[2*n];
            odd_r[n] = x_r[2*n+1];
            odd_i[n] = x_i[2*n+1];
        }
        FFT(even_r,even_i,even_FT_r,even_FT_i,N/2);
        FFT(odd_r,odd_i,odd_FT_r,odd_FT_i,N/2);
        
        for(k=0;k<N/2;++k)
        {
            w_r = cos(-k*2*M_PI/N);
            w_i = sin(-k*2*M_PI/N);
            // printf("N=%d, %f + %f i\n", N, w_r, w_i);  發現計算上的bug(非語言)印出來檢查
            y_r[k] = even_FT_r[k] + (w_r*odd_FT_r[k] - w_i*odd_FT_i[k]);
            y_i[k] = even_FT_i[k] + (w_r*odd_FT_i[k] + w_i*odd_FT_r[k]);
            y_r[k+N/2] = even_FT_r[k] - (w_r*odd_FT_r[k] - w_i*odd_FT_i[k]);
            y_i[k+N/2] = even_FT_i[k] - (w_r*odd_FT_i[k] + w_i*odd_FT_r[k]);

        }
        
        free(even_r);
        free(even_i);
        free(even_FT_r);
        free(even_FT_i);
        
        /*
        for(k=N/2;k<N-1;++k)
        {
            w_r = cos(-k*2*M_PI/N);
            w_i = sin(-k*2*M_PI/N);
            y_r[k] = even_r[k-N/2] + w_r*odd_r[k-N/2] - w_i*odd_i[k-N/2];
            y_i[k] = even_i[k-N/2] + w_r*odd_i[k-N/2] + w_i*odd_r[k-N/2];
        }
        
        u_r = (double *) malloc(N*sizeof(double));
        u_i = (double *) malloc(N*sizeof(double));
        v_r = (double *) malloc(N*sizeof(double));
        v_i = (double *) malloc(N*sizeof(double));
        for(n=0;n<N/2;n++)
        {
            u_r[n] = x_r[2*n];
            u_i[n] = x_i[2*n];
            u_r[n+N/2] = x_r[2*n+1];
            u_i[n+N/2] = x_i[2*n+1];
        }
        FFT(u_r, u_i, v_r,v_i, N/2);
        FFT(u_r+N/2, u_i+N/2, v_r+N/2, v_i+N/2, N/2);
         */
        return 0;
    }
    
    
