#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include<limits.h>
void recursion(long long int N,long long int a[][N],long long int i,long long int j,long long int b[],long long int val)
{

    if(i==j)
    {
        return;
    }
    long long int s=a[i][j];
    // printf("%d ",s);
    long long int left=i;
    long long int end=j;
    long long int mid=s;
    b[s]=val;
    // printf("1");
    recursion(N,a,left,mid-1,b,s);
    recursion(N,a,mid,end,b,s);
    
}
void merge(long long int *a, long long int *b,long long int *c,long long int *d, long long int l, long long int m, long long int r)
{

    long long int h, i, j, k;
    h = l;
    i = l;
    j = m + 1;
    while ((h <= m) && (j <= r))
    {

        if (a[h] <= a[j])
        {
            b[i] = a[h];
            d[i] = c[h];
            h++;
        }
        else
        {
            b[i] = a[j];
            d[i] = c[j];
            j++;
        }
        i++;
    }
    if (m < h)
    {
        for (k = j; k <= r; k++)
        {   
            b[i] = a[k];
            d[i] = c[k];
            i++;
        }
    }
    else
    {
        for (k = h; k <= m; k++)
        {
            b[i] = a[k];
            d[i] = c[k];
            i++;
        }
    }

    for (k = l; k <= r; k++)
    {
        a[k] = b[k];
        c[k] = d[k];
    }
}

void mergeSort(long long int *a, long long int *b,long long int *c,long long int *d ,long long int l, long long int r)
{

    long long int m;
    if (l < r)
    {

        m = (l + r) / 2;
        mergeSort(a, b,c,d ,l, m);
        mergeSort(a, b,c,d ,(m + 1), r);
        merge(a, b,c,d ,l, m, r);
    }
}
long long int sum(long long int freq[], long long int i, long long int j)
{
	long long int s = 0;
	for (long long int k = i; k <=j; k++)
	s += freq[k];
	return s;
}
long long int main(int argc, char *argv[])
{
    
    MPI_Init(&argc, &argv);
    int N;
    int size,rank;
    srand(time(NULL));
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Barrier( MPI_COMM_WORLD );
    double tbeg = MPI_Wtime();
    if (rank == 0)
    {
        int N;
        scanf("%d", &N);
        MPI_Bcast(&N, 1, MPI_LONG_LONG_INT,0, MPI_COMM_WORLD);
        long long int key[N];
        long long int freq[N];
        long long int parent[N][N];
        for(long long int i=0;i<N;i++)
        {
            for(long long int j=0;j<N;j++)
            {
               parent[i][j]=0;
               if(i==j)
               {
                  parent[i][j]=i;
               }
            }
           
        }
        for (long long int i = 0; i < N; i++)
        {
            scanf("%lld %lld", &key[i], &freq[i]);
        }
        // printf("\n");
        long long int length1 = (N) / size;
        long long int length2 = (N) - (size - 1) * length1;
        long long int* key1;
        long long int* freq1;
        freq1 = (long long int*)malloc(length1 * sizeof(long long int));
        key1 = (long long int*)malloc(length1 * sizeof(long long int));
        
        
        mergeSort(key, key1,freq,freq1, 0, length1 - 1);
        free(key1);
        free(freq1);
        if(size!=1)
        {

            for (long long int i = 1; i < size - 1; i++)
            {

                MPI_Send(&length1, 1, MPI_LONG_LONG_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(key + length1 * i, length1, MPI_LONG_LONG_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(freq + length1 * i, length1, MPI_LONG_LONG_INT, i, 0, MPI_COMM_WORLD);
            }
            MPI_Send(&length2, 1, MPI_LONG_LONG_INT, size - 1, 0, MPI_COMM_WORLD);
            MPI_Send(key + length1 * (size - 1), length2, MPI_LONG_LONG_INT, size - 1, 0, MPI_COMM_WORLD);
            MPI_Send(freq + length1 * (size - 1), length2, MPI_LONG_LONG_INT, size - 1, 0, MPI_COMM_WORLD);
            for (long long int i = 1; i < size - 1; i++)
            {
                MPI_Recv(key + length1 * i, length1, MPI_LONG_LONG_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(freq + length1 * i, length1, MPI_LONG_LONG_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Recv(key + length1 * (size - 1), length2, MPI_LONG_LONG_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(freq + length1 * (size - 1), length2, MPI_LONG_LONG_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


        }
        long long int* key2;
        long long int* freq2;
        freq2 = (long long int*)malloc(N * sizeof(long long int));
        key2 = (long long int*)malloc(N * sizeof(long long int));
        

        long long int len= N/size;
        if(len==0)
        {
            for(long long int i=0;i<N;i++)
            {
                key2[i]=key[i];
                freq2[i]=freq[i];
            }
        }
        else
        {
            long long int pointers[size];
            for(long long int i=0;i<size;i++)
            {
                pointers[i]=i*len;
                
            }
            for(long long int i=0;i<N;i++)
            {
                long long int mini=LLONG_MAX;
                long long int Pindex;
                long long int kindex;
                for(long long int j=0;j<size;j++)
                {
                    if(j!=size-1)
                    {
                        if(pointers[j]<(j+1)*len)
                        {
                            if(key[pointers[j]] < mini)
                            {
                                mini=key[pointers[j]];
                                Pindex=j;
                                kindex=pointers[j]; 
        
                            } 

                        }

                    }
                    else
                    {
                        if(pointers[j]<N)
                        {
                            if(key[pointers[j]] < mini)
                            {
                                mini=key[pointers[j]];
                                Pindex=j;
                                kindex=pointers[j];
                            } 

                        }
                        

                    }
                
                
                }
                key2[i]=key[kindex];
                freq2[i]=freq[kindex];
                pointers[Pindex]++;
            
            }

        }
        
        // free(key2);
        // free(freq2);
       
        long long int cost[N][N];
        for(long long int i=0;i<N;i++)
        {
            for(long long int j=0;j<N;j++)
            {
                cost[i][j]=0;
            }
        }
        for(long long int i=0;i<N;i++)
        {
            cost[i][i]=freq2[i];
        }
        for (long long int i = 1; i < size; i++)
        {
            MPI_Send(freq2,N, MPI_LONG_LONG_INT, i, 0, MPI_COMM_WORLD);
        }
        for (long long int L=2; L<=N; L++)
        {
	
            for (long long int i=0; i<N-L+1; i++)
            {
                
                // Get column number j from row number i and
                // chain length L
                long long int j = i+L-1;
                long long int off_set_sum = sum(freq2, i, j);
                // if(i==0 && j==2)
                // {
                //     printf("offset sum %d\n",off_set_sum);
                // }
                cost[i][j] = LLONG_MAX;
                long long int cout=(j-i+1)/size;
                long long int index=-1;
                long long int index1;

                if(cout==0)
                {
                    long long int start=i;
                    long long int end=j;
                    
                    long long int costs=LLONG_MAX;
                    for (long long int r=start; r<=end; r++)
                    {
                        // c = cost when keys[r] becomes root of this subtree
                        long long int c = ((r > i)? cost[i][r-1]:0) +
                                ((r < j)? cost[r+1][j]:0) +
                                off_set_sum;
                        
                        if (c < costs)
                        {
                            costs = c;
                            index=r;
                        } 

                    }
                    cost[i][j]=costs;
                    parent[i][j]=index;
                    // printf("%d ",index);
                    long long int reduction_result;
                    MPI_Allreduce(&costs, &reduction_result, 1, MPI_LONG_LONG_INT, MPI_MIN,MPI_COMM_WORLD);

                }
                else
                {
                    long long int start=i+rank*cout;
                    long long int end;
                    if(rank==size-1)
                    {
                        end=j;
                        
                    }
                    else
                    {
                        // end=i+((rank+1)*cout)-1;
                        end=start+cout-1;
                    }
                    long long int costs=LLONG_MAX;
                    for (long long int r=start; r<=end; r++)
                    {
                        // c = cost when keys[r] becomes root of this subtree
                        long long int c = ((r > i)? cost[i][r-1]:0) +
                                ((r < j)? cost[r+1][j]:0) +
                                off_set_sum;
                        if (c < costs)
                        {
                            costs = c;
                            index=r;
                        }
                            
                    }
                    long long int reduction_result;
                    MPI_Allreduce(&costs, &reduction_result, 1, MPI_LONG_LONG_INT, MPI_MIN,MPI_COMM_WORLD);
                    MPI_Barrier( MPI_COMM_WORLD );
                    // MPI_Barrier( MPI_COMM_WORLD );
                    cost[i][j]=reduction_result;
                    long long int spr=LLONG_MAX;
                    if(costs==reduction_result)
                    {
                        spr=index;
                    }
                    long long int red;
                    MPI_Allreduce(&spr, &red, 1, MPI_LONG_LONG_INT, MPI_MIN,MPI_COMM_WORLD);
                    parent[i][j]=red;
                   

                }
               
               
            }
          
            
        }
        
        long long int roots[N+1];
        for(long long int i=0;i<N+1;i++)
        {
            roots[i]=-1;
        }
        long long int dummy[N+1][N+1];
        for(long long int i=0;i<N+1;i++)
        {
            for(long long int j=0;j<N+1;j++)
            {
                dummy[i][j]=0;
                if(i==j) dummy[i][j]=i;
                
            }
            
        }
        for(long long int i=0;i<N;i++)
        {
            for(long long int j=0;j<N;j++)
            {
                
                dummy[i][j+1]=parent[i][j]+1;
  
                
            }
        }
        for(long long int i=0;i<N+1;i++)
        {
            dummy[i][i]=0;
        }
        recursion(N+1,dummy,0,N,roots,0);
        printf("%lld",cost[0][N-1]);
        printf("\n");
        for(long long int i=1;i<N+1;i++)
        {
            if(roots[i]==0)
            {
                long long int x=0;
                printf("%lld ",x);
            }
            else
            {
                printf("%lld ",key2[roots[i]-1]);
            }
            
        }
        printf("\n");
        free(key2);
        free(freq2);


    }

    if(rank!=0)
    {
        long long int length2;
        MPI_Recv(&length2, 1, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Bcast(&N, 1, MPI_LONG_LONG_INT,0, MPI_COMM_WORLD);
        long long int key[length2];
        long long int freq[length2];
        long long int* key1;
        long long int* freq1;
        freq1 = (long long int*)malloc(length2 * sizeof(long long int));
        key1 = (long long int*)malloc(length2 * sizeof(long long int));
        MPI_Recv(key, length2, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(freq, length2, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        mergeSort(key, key1,freq,freq1, 0, length2 - 1);
        free(key1);
        free(freq1);
        MPI_Send(key, length2, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(freq, length2, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);

        long long int fre[N];

        MPI_Recv(fre,N, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        long long int *data = (long long int *)malloc(N*N*sizeof(long long int));
        long long int **cost= (long long int **)malloc(N*sizeof(long long int*));
        for (long long int i=0; i<N; i++)
            cost[i] = &(data[N*i]);

        for(long long int i=0;i<N;i++)
        {
            for(long long int j=0;j<N;j++)
            {
                cost[i][j]=0;
            }
        }
        for(long long int i=0;i<N;i++)
        {
            cost[i][i]=fre[i];
           
        }
        
        for (long long int L=2; L<=N; L++)
        {
	
            for (long long int i=0; i<N-L+1; i++)
            {
                long long int j = i+L-1;
                long long int off_set_sum = sum(fre, i, j);
                cost[i][j] = LLONG_MAX;
                long long int cout=(j-i+1)/size;
                if(cout!=0)
                {
                    long long int start=i+rank*cout;
                    long long int end;
                    if(rank==size-1)
                    {
                        end=j;
                        
                    }
                    else
                    {
                        
                        end=start+cout-1;
                    }
                    long long int costs=LLONG_MAX;
                    long long int index=-1;
                    for (long long int r=start; r<=end; r++)
                    {
                        
                        long long int c = ((r > i)? cost[i][r-1]:0) +
                                ((r < j)? cost[r+1][j]:0) +
                                off_set_sum;
                        
                
                        if (c < costs)
                        {
                            costs = c;
                            index = r;
    
                        }
                            
                    }
                    long long int reduction_result;
                    MPI_Allreduce(&costs, &reduction_result, 1, MPI_LONG_LONG_INT, MPI_MIN,MPI_COMM_WORLD);
                    MPI_Barrier( MPI_COMM_WORLD );
                    long long int spr=LLONG_MAX;
                    if(costs==reduction_result)
                    {
                        spr=index;
                    }
                    long long int red;
                    MPI_Allreduce(&spr, &red, 1, MPI_LONG_LONG_INT, MPI_MIN,MPI_COMM_WORLD);
                   
                    cost[i][j]=reduction_result;
        

                }
                else
                {
                    long long int costs=LLONG_MAX;
                    long long int reduction_result;
                    MPI_Allreduce(&costs, &reduction_result, 1, MPI_LONG_LONG_INT, MPI_MIN,MPI_COMM_WORLD);
                    cost[i][j]=reduction_result;
                }
                
            }
        }
    }



    MPI_Barrier( MPI_COMM_WORLD );
    double elapsedTime = MPI_Wtime() - tbeg;
    double maxTime;
    MPI_Reduce( &elapsedTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    if ( rank == 0 ) {
        printf( "Total time (s): %f\n", maxTime );
    }  
    MPI_Finalize();

}