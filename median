//
//  main.m
//  sort.v.1
//
//  Created by Nancy Wood on 2015/4/19.
//  Copyright (c) 2015年 Nancy Wood. All rights reserved.
//

#import <Foundation/Foundation.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, const char * argv[]) {
    @autoreleasepool {
        // insert code here...
        int i, n, *v, temp, j;
        srand(time(0));
        n = 10;
        v = (int *) malloc( n * sizeof(int));
        for(i=0; i<n; ++i)
        {
            v[i] = rand() % 100;
            printf("%d,",v[i]);
        }
        printf("\n");
        
        
        // 選擇排序
        // 選最大值，放在位置n
        // 剩下n-1個選最大值，放在位置n-1
        // 到剩下1個停
        
        for(i=0; i<n; ++i)
        {
            int max=i;
            
            
                for(j=i; j<n; ++j)
                {
                    if(v[j]>v[max])
                        max=j;
                }
                if(max!=i)
                {
                    temp = v[max];
                    v[max]=v[i];
                    v[i]=temp;
                }
                for(int k=0; k<n; ++k)
                {
                    printf("%d,",v[k]);
                }
                printf("\n");
            
        }
        
        // median
        // 因為n已經知道是10，所以median是第5、6大的數的平均
        // 是說我們也不會對未知的n做排序吼...因為會排不完= ="
        // 可以玩有上界的random n個數排序，就可以判斷奇偶數再決定median算法嗎
        // 嗚嗚等急迫的事忙完再回來想
        
        int a=v[4], b=v[5];
        int sum=a+b;
        double k=sum/2.0;  //2.0太神奇了啦XDDD
        printf("median is %f",k);
        printf("\n");
        
        
        
        
        //NSLog(@"Hello, World!");
    }
    return 0;
}

