/*

The purpose of this file is to accelerate the execution of the ribbon carving algorithm.

The Mex function in this file passes the arguments from the MATLAB function: N2=ribboncarvemainC(X,C,flex,eps,vRibbonCount,hRibbonCount);

Compile this file in MATLAB to generate a MEX-file before running.

Variables:
	input:
		X: video data (4D data (y,x,rgb,t) using X_INDEX() for indexing) (same variable as output)
		C: cost data (3D data (y,x,t)) (same variable as output)
		flex: flex parameter (integer, >= 1)
		eps: maximum allowed cost of removing a ribbon for stopping criterion
		vRibbonCount: the number of vertical ribbons removed
		hRibbonCount: the number of horizontal ribbons removed
	output:
		N2: the remaining number of frames after ribbon carving

Huan-Yu Wu

*/


#include "mex.h"

//index conversion from multi-dimension to 1-D array using #define
#define TWO_D_INDEX(y,x,height) ((y)+(height)*(x))
#define THREE_D_INDEX(y,x,t,height,weight) ((y)+(height)*(x)+(height)*(weight)*(t))
#define X_INDEX(y,x,rgb,t,height,weight) ((y)+(height)*(x)+(height)*(weight)*(rgb)+(height)*(weight)*3*(t))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	unsigned char *X, *C;
	double *flex, *eps, *vRibbonCount, *hRibbonCount, *N2;
	int *size, Ch, Cv, H, W, N, minIndv, minIndh;
    unsigned short int *Mx, *My;  //2D matrices storing cumulative minimum costs in horizontal and vertical projections
	char *pathMx, *pathMy;  //2D matrices storing directions of paths in horizontal and vertical projections

    
	void calculatePathWithMinCost(unsigned short int *, char *, unsigned char *, double *, char, int, int, int);
	void carveRibbon(unsigned char *, unsigned char *, char *, int, char, int, int, int);
	int findMin(unsigned short int *, int, int, int, int *, int);

	//Input variables
	X = mxGetData(prhs[0]);
	C = mxGetData(prhs[1]);
	size = mxGetDimensions(prhs[1]);
	flex = mxGetPr(prhs[2]);
	eps = mxGetPr(prhs[3]);
	vRibbonCount = mxGetPr(prhs[4]);
	hRibbonCount = mxGetPr(prhs[5]);

	//Calculations
	H = size[0];
    W = size[1];
    N = size[2];  //number of frames

	if (!(My = (unsigned short int *)malloc(H * N * sizeof(unsigned short int))))
		printf("Memory allocation failed for My in MEX.\n");
	if (!(pathMy = (char *)malloc(H * N * sizeof(char))))
		printf("Memory allocation failed for pathMy in MEX.\n");
	if (!(Mx = (unsigned short int *)malloc(W * N * sizeof(unsigned short int))))
		printf("Memory allocation failed for Mx in MEX.\n");
	if (!(pathMx = (char *)malloc(W * N * sizeof(char))))
        printf("Memory allocation failed for pathMx in MEX.\n");

	calculatePathWithMinCost(My, pathMy, C, flex, 'v', H, W, N);
	calculatePathWithMinCost(Mx, pathMx, C, flex, 'h', H, W, N);

    //Cv/Ch: minimum cost of removing a vertical/horizontal ribbon
	Cv = findMin(My, H-1, 0, N-1, &minIndv, H);
	Ch = findMin(Mx, W-1, 0, N-1, &minIndh, W);

	while (Ch <= *eps || Cv <= *eps)
	{
		if (Cv < Ch)			//remove a vertical ribbon
		{
			carveRibbon(X, C, pathMy, minIndv, 'v', H, W, N--);
			(*vRibbonCount)++;
		}
		else if (Ch < Cv)		//remove a horizontal ribbon
		{
			carveRibbon(X, C, pathMx, minIndh, 'h', H, W, N--);
			(*hRibbonCount)++;
		}
		else				//remove a vertical ribbon or a horizontal ribbon (choose by random)
		{
			if (rand() % 2)
			{
				carveRibbon(X, C, pathMy, minIndv, 'v', H, W, N--);
				(*vRibbonCount)++;
			}
			else
			{
				carveRibbon(X, C, pathMx, minIndh, 'h', H, W, N--);
				(*hRibbonCount)++;
			}
		}

		//recalculate the cumulative minimum cost
		calculatePathWithMinCost(My, pathMy, C, flex, 'v', H, W, N);
		calculatePathWithMinCost(Mx, pathMx, C, flex, 'h', H, W, N);
		Cv = findMin(My, H-1, 0, N-1, &minIndv, H);
		Ch = findMin(Mx, W-1, 0, N-1, &minIndh, W);
        
	}
	free(My);
	free(pathMy);
	free(Mx);
	free(pathMx);

	//Output variable
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	N2 = mxGetPr(plhs[0]);
	*N2 = (double)N;

}

void calculatePathWithMinCost(unsigned short int *M, char *pathM, unsigned char *C, double *flex, char direction, int H, int W, int N)
{
    // In this function, we project cost buffer on either vertical or horizontal direction, then
    // use dynamic programming to calculate path with minimum cost in each direction.
    //
    //   M: intermediate 2D matrix for storing cumulative cost
    // Input data:
    //   C: 3D cost buffer
    //   flex: flex parameter
    //   direction: 'v' for vertical; 'h' for horizontal
    //   H: height of buffer
    //   W: width of buffer
    //   N: number of frames in the buffer
    // Output data:
    //   pathM: matrix that will be used later to trace back to top row.
    //          0: same column; -1: one column left; 1: one column right; -2: two columns left; etc.
    
	register unsigned int indexM, indexC;
	unsigned short int tempM;
	int i, j, k, minInd;
    
	switch (direction)
	{
	case 'v':

		//calculate the values in the first row of M (size of M: H*N)
		for (j = 0; j < N; ++j)
		{
			indexM = TWO_D_INDEX(0, j, H);
			indexC = H * W * j;

            //projection
			tempM = 0;
			for (k = 0; k < W; ++k)
				tempM += C[indexC + H * k];

			M[indexM] = tempM;
		}

		//calculate the values from the second row to the last row of M
		for (i = 1; i < H; ++i)
		{
			for (j = 0; j < N; ++j)
			{
				indexM = TWO_D_INDEX(i, j, H);
				indexC = i + H * W * j;

                //projection
				tempM = 0;
				for (k = 0; k < W; ++k)
					tempM += C[indexC + H * k];

                //add minimum cumulative cost from upper row and record the path
				if (j < *flex)
				{
					tempM += findMin(M, i - 1, 0, j + (int)(*flex), &minInd, H);
					pathM[indexM] = minInd - j;
				}
				else if (j < N - *flex)
				{
					tempM += findMin(M, i - 1, j - (int)(*flex), j + (int)(*flex), &minInd, H);
					pathM[indexM] = minInd - *flex;
				}
				else
				{
					tempM += findMin(M, i - 1, j - (int)(*flex), N - 1, &minInd, H);
					pathM[indexM] = minInd - *flex;
				}

				M[indexM] = tempM;
			}
		}

		break;

	case 'h':

		//calculate the values in the first row of M (size of M: W*N)
		for (j = 0; j < N; ++j)
		{
			indexM = TWO_D_INDEX(0, j, W);
			indexC = H * W * j;

            //projection
			tempM = 0;
			for (k = 0; k < H; ++k)
				tempM += C[indexC + k];

			M[indexM] = tempM;
		}

		//calculate the values from the second row to the last row of M
		for (i = 1; i < W; ++i)
		{
			for (j = 0; j < N; ++j)
			{
				indexM = TWO_D_INDEX(i, j, W);
				indexC = H * i + H * W * j;

                //projection
				tempM = 0;
				for (k = 0; k < H; ++k)
					tempM += C[indexC + k];

                //add minimum cumulative cost from upper row and record the path
				if (j < *flex)
				{
					tempM += findMin(M, i - 1, 0, j + (int)(*flex), &minInd, W);
					pathM[indexM] = minInd - j;
				}
				else if (j < N - *flex)
				{
					tempM += findMin(M, i - 1, j - (int)(*flex), j + (int)(*flex), &minInd, W);
					pathM[indexM] = minInd - *flex;
				}
				else
				{
					tempM += findMin(M, i - 1, j - (int)(*flex), N - 1, &minInd, W);
					pathM[indexM] = minInd - *flex;
				}

				M[indexM] = tempM;
			}
		}

		break;


	/* //ORIGINAL MATLAB METHOD
	case 'v':
		//sum up the cost along the W direction, and save it to M (size of M: H*N)
		for(i=0;i<H;i++)
		{
			for(j=0;j<N;j++)
			{
				M[TWO_D_INDEX(i,j,H)] = 0;
				for(k=0;k<W;k++)
					M[TWO_D_INDEX(i,j,H)] += C[THREE_D_INDEX(i,k,j,H,W)];
			}
		}

		//compute the cumulative cost
		for(i=1;i<H;i++)
		{
			for(j=*flex;j<N-*flex;j++)
			{
				M[TWO_D_INDEX(i,j,H)] += min(M, i-1, j-(int)(*flex), j+(int)(*flex), minInd, H);
				pathM[TWO_D_INDEX(i,j,H)] = *minInd - *flex;
				//pathM[TWO_D_INDEX(i,j,H)] = *minInd;
			}
			for(j=0;j<*flex;j++)
			{
				M[TWO_D_INDEX(i,j,H)] += min(M, i-1, 0, j+(int)(*flex), minInd, H);
				pathM[TWO_D_INDEX(i,j,H)] = *minInd - j;
				//pathM[TWO_D_INDEX(i,j,H)] = *minInd;
			}
			for(j=N-*flex;j<N;j++)
			{
				M[TWO_D_INDEX(i,j,H)] += min(M, i-1, j-(int)(*flex), N-1, minInd, H);
				pathM[TWO_D_INDEX(i,j,H)] = *minInd - *flex;
				//pathM[TWO_D_INDEX(i,j,H)] = *minInd;
			}
		}

		break;

	case 'h':
		//sum up the cost along the H direction, and save it to M (size of M: W*N)
		for(i=0;i<W;i++)
		{
			for(j=0;j<N;j++)
			{
				M[TWO_D_INDEX(i,j,W)] = 0;
				for(k=0;k<H;k++)
					M[TWO_D_INDEX(i,j,W)] += C[THREE_D_INDEX(k,i,j,H,W)];
			}
		}

		//compute the cumulative cost
		for(i=1;i<W;i++)
		{
			for(j=*flex;j<N-*flex;j++)
			{
				M[TWO_D_INDEX(i,j,W)] += min(M, i-1, j-(int)(*flex), j+(int)(*flex), minInd, W);
				pathM[TWO_D_INDEX(i,j,W)] = *minInd - *flex;
			}
			for(j=0;j<*flex;j++)
			{
				M[TWO_D_INDEX(i,j,W)] += min(M, i-1, 0, j+(int)(*flex), minInd, W);
				pathM[TWO_D_INDEX(i,j,W)] = *minInd - j;
			}
			for(j=N-*flex;j<N;j++)
			{
				M[TWO_D_INDEX(i,j,W)] += min(M, i-1, j-(int)(*flex), N-1, minInd, W);
				pathM[TWO_D_INDEX(i,j,W)] = *minInd - *flex;
			}
		}

		break;
		*/
	}
}


void carveRibbon(unsigned char *X, unsigned char *C, char *pathM, int p, char direction, int H, int W, int N)
{
    // In this function, We remove a ribbon from video and cost buffer.
    // After this function, the number of frames will be reduced by one.
    //
    // Input:
    //   X: 4D video buffer (same variable as output)
    //   C: 3D cost buffer (same variable as output)
    //   pathM: 2D matrix that stores directions to trace back
    //   p: current column location of path (starting from last row)
    //   direction: 'v' for vertical; 'h' for horizontal
    //   H: height of buffer
    //   W: width of buffer
    //   N: number of frames in the buffer
    
    
	int i, j, k;
	register unsigned int indexX, indexC, HW=H*W;

	switch (direction)
	{
	case 'v':

		for (i = H - 1; i >= 0; --i)
		{
			for (j = 0; j < W; ++j)
			{
				indexX = X_INDEX(i, j, 0, 0, H, W) + HW * 3 * p;			//temporary index value for X
				indexC = THREE_D_INDEX(i, j, 0, H, W) + HW * p;		//temporary index value for C
				for (k = p; k < N - 1; ++k)
				{
					X[indexX] = X[indexX + HW * 3];					//X[X_INDEX(i,j,0,k,H,W)]=X[X_INDEX(i,j,0,k+1,H,W)];
					X[indexX + HW] = X[indexX + HW * 4];			//X[X_INDEX(i,j,1,k,H,W)]=X[X_INDEX(i,j,1,k+1,H,W)];
					X[indexX + HW * 2] = X[indexX + HW * 5];			//X[X_INDEX(i,j,2,k,H,W)]=X[X_INDEX(i,j,2,k+1,H,W)];
					C[indexC] = C[indexC + HW];					//C[THREE_D_INDEX(i,j,k,H,W)]=C[THREE_D_INDEX(i,j,k+1,H,W)];
					indexX += HW * 3;
					indexC += HW;

				}
			}
			p += pathM[TWO_D_INDEX(i, p, H)];			//calculate ribbon path
		}

		break;

	case 'h':

		for (i = W - 1; i >= 0; --i)
		{
			for ( j = 0; j < H; ++j)
			{
				indexX = X_INDEX(j, i, 0, 0, H, W) + HW * 3 * p;
				indexC = THREE_D_INDEX(j, i, 0, H, W) + HW * p;
				for ( k = p; k < N - 1; ++k)
				{
					X[indexX] = X[indexX + HW * 3];					//X[X_INDEX(j,i,0,k,H,W)]=X[X_INDEX(j,i,0,k+1,H,W)];
					X[indexX + HW] = X[indexX + HW * 4];			//X[X_INDEX(j,i,1,k,H,W)]=X[X_INDEX(j,i,1,k+1,H,W)];
					X[indexX + HW * 2] = X[indexX + HW * 5];			//X[X_INDEX(j,i,2,k,H,W)]=X[X_INDEX(j,i,2,k+1,H,W)];
					C[indexC] = C[indexC + HW];					//C[THREE_D_INDEX(j,i,k,H,W)]=C[THREE_D_INDEX(j,i,k+1,H,W)];
					indexX += HW * 3;
					indexC += HW;

				}
			}
			p += pathM[TWO_D_INDEX(i, p, W)];
		}

		break;

	}
}




int findMin(unsigned short int *M, int row, int startCol, int endCol, int *minInd, int h)
{
    // This function finds the minimum value within a range in a row of a given matrix.
    //
    // Input:
    //   M: 2D matrix of cumulative cost
    //   row: row in the matrix
    //   startCol: starting column
    //   endCol: ending column
    // Output:
    //   minInd: index of minimum value (startCol being 0)
    //   tempMin: minimum value
    
	register int i;
	int tempMin;
	
	tempMin = M[TWO_D_INDEX(row, startCol, h)];
	*minInd = 0;
	for (i = startCol + 1; i <= endCol; ++i)
	{
		if (M[TWO_D_INDEX(row, i, h)] < tempMin)
		{
			tempMin = M[TWO_D_INDEX(row, i, h)];
			*minInd = i - startCol;
		}
	}

	return tempMin;

}

