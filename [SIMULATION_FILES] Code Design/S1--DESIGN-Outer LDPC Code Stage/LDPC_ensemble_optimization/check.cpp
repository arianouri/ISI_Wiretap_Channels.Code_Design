#include <mex.h>
#include <math.h>

void Check( double  *pa, double *pb, double *pc, unsigned int *look_up, unsigned int N);

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
	double *pa;
	double *pb;
	double *pc;
	unsigned int *look_up;
    unsigned int N;

        pa = mxGetPr(prhs[0]);
	pb = mxGetPr(prhs[1]);
        look_up = (unsigned int *)mxGetData(prhs[2]);
        N = (unsigned int)mxGetScalar(prhs[3]);
	plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
	pc = mxGetPr(plhs[0]);

	Check(pa,pb,pc,look_up,N);
	
}


void Check( double  *pa, double *pb, double *pc, unsigned int *look_up, unsigned int N)
{
   unsigned int i,j,k;
   unsigned int c=0;
   unsigned int t=(N-1)/2;
   for (i=0;i<N;i++)
   {
        pc[i]=0.0;
   }

   pc[t]=pa[t]+pb[t]-pa[t]*pb[t];

   for (j=1;j<=t;j++)
	{
        for (i=j;i<=t;i++)
                {
                k=look_up[c];
                c++;
                if (i==j)
                        {
                        pc[k]+=pa[t+j]*pb[t+i];
                        pc[k]+=pa[t-j]*pb[t-i];
                        pc[N-1-k]+=pa[t+j]*pb[t-i];
                        pc[N-1-k]+=pa[t-j]*pb[t+i];
                        }
                else
                        {
                        pc[k]+=pa[t+j]*pb[t+i];
                        pc[k]+=pb[t+j]*pa[t+i];
                        pc[k]+=pa[t-j]*pb[t-i];
                        pc[k]+=pb[t-j]*pa[t-i];
                        pc[N-1-k]+=pa[t+j]*pb[t-i];
                        pc[N-1-k]+=pb[t+j]*pa[t-i];
                        pc[N-1-k]+=pa[t-j]*pb[t+i];
                        pc[N-1-k]+=pb[t-j]*pa[t+i];
                        }
                }
        }



}


