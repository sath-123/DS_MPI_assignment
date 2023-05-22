#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <complex.h>

void verify(float _Complex arr[],int answer[],int length,int K)
{
    
    for(int j=0;j<length;j++)
    {
        int i=0;
        float Im=0;
        float Re=0;
        while (((i < K+1) && (Im*Im) + (Re * Re) <4))
        {
            float temp = (Re * Re) - (Im * Im) + crealf(arr[j]);
            Im = 2.0 * Re * Im + cimagf(arr[j]);
            Re = temp;
            i++;
        }
        if(i>K)
        {
            answer[j]=1;
        }
        else
        {
            answer[j]=0;
        }
        
    }
}
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int size, rank,K;
    srand(time(NULL));
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);
    double tbeg = MPI_Wtime();
    
    if(rank==0)
    {
        int N,M;
        scanf("%d %d %d",&N,&M,&K);
        MPI_Bcast(&K, 1, MPI_INT,0, MPI_COMM_WORLD);
        float CRmax= 1;
        float CRmin= -1.5;
        float CImax= 1;
        float CImin= -1;
        int xl=N-1;
        int yl=M-1;
        float xlen=(CRmax-CRmin)/(float)xl;
        float ylen=(CImax-CImin)/(float)yl;
        float _Complex* complexToSend;
        complexToSend = (float _Complex*)malloc(N*M*sizeof(float _Complex));
        int* Answer;
        Answer = (int*)malloc(N*M*sizeof(int*));
        for(int i=0;i<N*M;i++)
        {
            int row=i%N;
            int col=i/M;
            float CX = CRmin+row*xlen;
            float CY=CImax-col*ylen;
            complexToSend[i] = CX+CY*I;
        }
        int length1=(N*M)/size;
        int length2=(N*M)-(size-1)*length1;

        int len=length1;
        verify(complexToSend,Answer,len,K);
        if(size!=1)
        {
            for(int i=1;i<size-1;i++)
            {
                MPI_Send(&len,1,MPI_INT,i,0,MPI_COMM_WORLD);
                MPI_Send(complexToSend+len*i, length1, MPI_C_COMPLEX,i, 0, MPI_COMM_WORLD);
            }
            MPI_Send(&length2,1,MPI_INT,size-1,0,MPI_COMM_WORLD);
            MPI_Send(complexToSend+len*(size-1), length2, MPI_C_COMPLEX,size-1, 0, MPI_COMM_WORLD);
            for(int i=1;i<size-1;i++)
            {
                MPI_Recv(Answer+len*i, length1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Recv(Answer+len*(size-1), length2, MPI_INT, size-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }
        
        int row=0;
        for(int i=0;i<N*M;i++)
        {
            if(row!=0 && i%N==0)
            {
                printf("\n");
            }
            if(i%N==0)
            {
                row =1;
            }
            printf("%d ",Answer[i]);
        }
        printf("\n");
        free(Answer);
        free(complexToSend);
        

    }
    if(rank!=0)
    {
        MPI_Bcast(&K, 1, MPI_INT,0, MPI_COMM_WORLD);
        int length2;
        MPI_Recv(&length2,1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        float _Complex* complexToSend;
        complexToSend = (float _Complex*)malloc(length2*sizeof(float _Complex));
        int* answer;
        answer=(int*)malloc(length2*sizeof(float _Complex));
        MPI_Recv(complexToSend,length2,MPI_C_COMPLEX, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int var=1;
        int va=var;
        verify(complexToSend,answer,length2,K);
        MPI_Send(answer, length2, MPI_INT,0, 0, MPI_COMM_WORLD);
        free(answer);
        free(complexToSend);
        
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


