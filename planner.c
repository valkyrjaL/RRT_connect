/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"

#include <float.h>

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]


/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

//the length of each link in the arm (should be the same as the one used in runtest.m)
#define LINKLENGTH_CELLS 10

typedef struct {
  int X1, Y1;
  int X2, Y2;
  int Increment;
  int UsingYIndex;
  int DeltaX, DeltaY;
  int DTerm;
  int IncrE, IncrNE;
  int XIndex, YIndex;
  int Flipped;
} bresenham_param_t;


void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size)
{
    double cellsize = 1.0;
	//take the nearest cell
	*pX = (int)(x/(double)(cellsize));
	if( x < 0) *pX = 0;
	if( *pX >= x_size) *pX = x_size-1;

	*pY = (int)(y/(double)(cellsize));
	if( y < 0) *pY = 0;
	if( *pY >= y_size) *pY = y_size-1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params)
{
  params->UsingYIndex = 0;

  if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
    (params->UsingYIndex)++;

  if (params->UsingYIndex)
    {
      params->Y1=p1x;
      params->X1=p1y;
      params->Y2=p2x;
      params->X2=p2y;
    }
  else
    {
      params->X1=p1x;
      params->Y1=p1y;
      params->X2=p2x;
      params->Y2=p2y;
    }

   if ((p2x - p1x) * (p2y - p1y) < 0)
    {
      params->Flipped = 1;
      params->Y1 = -params->Y1;
      params->Y2 = -params->Y2;
    }
  else
    params->Flipped = 0;

  if (params->X2 > params->X1)
    params->Increment = 1;
  else
    params->Increment = -1;

  params->DeltaX=params->X2-params->X1;
  params->DeltaY=params->Y2-params->Y1;

  params->IncrE=2*params->DeltaY*params->Increment;
  params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
  params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

  params->XIndex = params->X1;
  params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y)
{
  if (params->UsingYIndex)
    {
      *y = params->XIndex;
      *x = params->YIndex;
      if (params->Flipped)
        *x = -*x;
    }
  else
    {
      *x = params->XIndex;
      *y = params->YIndex;
      if (params->Flipped)
        *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params)
{
  if (params->XIndex == params->X2)
    {
      return 0;
    }
  params->XIndex += params->Increment;
  if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
    params->DTerm += params->IncrE;
  else
    {
      params->DTerm += params->IncrNE;
      params->YIndex += params->Increment;
    }
  return 1;
}



int IsValidLineSegment(double x0, double y0, double x1, double y1, double*	map,
		   int x_size,
 		   int y_size)

{
	bresenham_param_t params;
	int nX, nY; 
    short unsigned int nX0, nY0, nX1, nY1;

    //printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
    
	//make sure the line segment is inside the environment
	if(x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

    //printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	//iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do {
		get_current_point(&params, &nX, &nY);
		if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
            return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double*	map,
		   int x_size, int y_size)
{
    double x0,y0,x1,y1;
    int i;
    
 	//iterate through all the links starting with the base
	x1 = ((double)x_size)/2.0;
    y1 = 0;
	for(i = 0; i < numofDOFs; i++)
	{
		//compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
		y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

		//check the validity of the corresponding line segment
		if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
				return 0;
	}    
	return 1;
}

/******************************************************************/
typedef struct{
	struct Node_t* next;
	struct Node_t* parent;
	double* config;
} Node_t;

void print_config(double* config){
	for(int i=0; i<5; i++){
		printf("%f\t", config[i]);
	}printf("\n");
}

double** init2Darray(int C, int R){
	double** A = (double**) malloc(R*sizeof(double*));
	for (int i = 0; i < R; i++){
        A[i] = (double*) malloc(C*sizeof(double)); 
	}
	return A;
}

void RANDOM_CONFIG(double* q, int DOF){
	for(int i = 0; i < DOF; i++){
		q[i] = ((float)rand()/(float)(RAND_MAX)) *2*PI;
	}
}

double dist(double *config1, double *config2, int dof) {
	double distance = 0;
	for (int i = 0; i < dof; i++){
		if(distance < fabs(config1[i] - config2[i])){
			distance = fabs(config1[i] - config2[i]);
		}
    }
	return distance;
}

void findXYcoord(double* angles, double* rx, double* ry, int x_size, int numofDOFs){
	double x0,y0,x1,y1;
    int i;
    
 	//iterate through all the links starting with the base
	x1 = ((double)x_size)/2.0;
    y1 = 0;
	for(i = 0; i < numofDOFs; i++)
	{
		//compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
		y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);
	}
	*rx = x1;
	*ry = y1;
}

Node_t* NEAREST_NEIGHBOR(Node_t* Tree, double* qrand, int x_size, int dof){
	Node_t* nearNode = NULL;
	double distance = 0.0;
	double min = FLT_MAX;
	double xrand, yrand;
	double xnear, ynear;
	findXYcoord(qrand, &xrand, &yrand, x_size, dof);
	while(Tree != NULL){
		findXYcoord(Tree->config, &xnear, &ynear, x_size, dof);
		distance = (xnear-xrand)*(xnear-xrand) + (ynear-yrand)*(ynear-yrand);
		if(distance < min){
			nearNode = Tree;
			min = distance;
		}
		Tree = Tree->next;
	}
	return nearNode;
}

// EXTEND from qnear to qrand
int EXTEND(Node_t* Tree, double* qrand, double* qnew, int x_size, int y_size, int numofDOFs, double* map){
	Node_t* qnear = NEAREST_NEIGHBOR(Tree, qrand, x_size, numofDOFs);
	
	// printf("rand:\t");	print_config(qrand);
	// printf("near:\t");	print_config(qnear->config);
	Node_t* lastNode = Tree;
	while(lastNode->next != NULL){
		lastNode = lastNode->next;
	}
	
	double distance = 0;
	int i,j;
	for (j = 0; j < numofDOFs; j++){
		if(distance < fabs(qrand[j] - qnear->config[j]))
			distance = fabs(qrand[j] - qnear->config[j]);
	}
	int numofsamples = (int)(distance/(PI/20));
	if(numofsamples < 2){
		return 2;
	}
	for (i = 0; i < numofsamples; i++){
		for(j = 0; j < numofDOFs; j++){
			qnew[j] = qnear->config[j] + ((double)(i)/(numofsamples-1))*(qrand[j] - qnear->config[j]);
		}
		if(!IsValidArmConfiguration(qnew, numofDOFs, map, x_size, y_size)){
			if(i<2){
				// printf("Trapped\n");
				return 2;	// Trapped
			}else{
				// get back to previous state
				for(j = 0; j < numofDOFs; j++){
					qnew[j] = qnear->config[j] + ((double)(i-1)/(numofsamples-1))*(qrand[j] - qnear->config[j]);
				}
				// printf("Advanced\n");	// add to tree
				double* qaddv = malloc(numofDOFs * sizeof(double));
				memcpy(qaddv, qnew, numofDOFs * sizeof(double));
				Node_t* NodeN = malloc(sizeof(Node_t));
				NodeN->next = NULL;
				NodeN->parent = qnear;
				NodeN->config = qaddv;
				lastNode->next = NodeN;
				// printf("new a:\t");	print_config(qnew);
				return 1;		// Advanced
			}
		}
	}
	// printf("Reached\n");	// add to tree
	double* qreach = malloc(numofDOFs * sizeof(double));
	// memcpy(qreach, qrand, numofDOFs * sizeof(double));
	memcpy(qreach, qnew, numofDOFs * sizeof(double));
	Node_t* NodeN = malloc(sizeof(Node_t));
	NodeN->next = NULL;
	NodeN->parent = qnear;
	NodeN->config = qreach;
	lastNode->next = NodeN;
	// printf("new r:\t");	print_config(qnew);
	return 0;	// Reached
}

void SWAP(Node_t** A, Node_t** B){
	Node_t* C = *A;
	*A = *B;
	*B = C;
}

// return the last node and set the length of tree at N
Node_t* checkTree(Node_t* Tree, int* N){
	int i = 0;
	while(Tree->next != NULL){
		i++;
		Tree = Tree->next;
	}
	
	*N = i;
	return Tree;
}

Node_t* checkConfig(Node_t* Tree, double* config_where){
	while(Tree->next != NULL){
		if(Tree->config == config_where){
			printf("Config is in this tree");
			break;
		}
		Tree = Tree->next;
	}
	
	return Tree;
}

int trackTree(Node_t* Tree){
	int i = 1;
	while(Tree->parent != NULL){
		i++;
		Tree = Tree->parent;
	}
	return i;
}

static void planner(
		   double*	map,
		   int x_size,
 		   int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
	   int numofDOFs,
	   double*** plan,
	   int* planlength)
{
	//no plan by default
	*plan = NULL;
	*planlength = 0;
	
    const int K = 10000;
	Node_t* TreeA = malloc(sizeof(Node_t));
	TreeA->next = NULL;
	TreeA->parent = NULL;
	TreeA->config = armstart_anglesV_rad;
	Node_t* TreeB = malloc(sizeof(Node_t));
	TreeB->next = NULL;
	TreeB->parent = NULL;
	TreeB->config = armgoal_anglesV_rad;

	double* qrand = malloc(numofDOFs * sizeof(double));
	double* qnew = malloc(numofDOFs * sizeof(double));
	
	for(int k = 1; k<K; k++){
		RANDOM_CONFIG(qrand, numofDOFs);
		if(IsValidArmConfiguration(qrand, numofDOFs, map, x_size, y_size)){
			if(!(EXTEND(TreeA, qrand, qnew, x_size, y_size, numofDOFs, map) == 2)){
				memcpy(qrand, qnew, numofDOFs * sizeof(double));
				if(EXTEND(TreeB, qrand, qnew, x_size, y_size, numofDOFs, map) == 0){	//0	Reached
					printf("Found gaol!!!!!!!!! K = %d \n", k);
					//********Backtrack for path********//
					int NA, NB;
					Node_t* LastA = checkTree(TreeA, &NA);
					Node_t* LastB = checkTree(TreeB, &NB);
					printf("TreeA: %d\t", NA);
					printf("TreeB: %d\n", NB);
					NA = trackTree(LastA);
					NB = trackTree(LastB);
					printf("A path: %d\t", NA);
					printf("B path: %d\n", NB);
					
					int N  = NA + NB;
					int NC = 0;
					*plan = init2Darray(numofDOFs, N);
					Node_t* TreeC = NULL;
					if(TreeA->config == armstart_anglesV_rad){
						TreeC = LastA;
						NC = NA;
					}else{
						TreeC = LastB;
						NC = NB;
					}
					for(int k=0; k<NC; k++){
						memcpy((*plan)[NC-1-k], TreeC->config, numofDOFs * sizeof(double));
						TreeC = TreeC -> parent;
					}
					if(TreeA->config == armstart_anglesV_rad){TreeC = LastB;}else{TreeC=LastA;}
					for(int k=NC; k<N; k++){
						memcpy((*plan)[k], TreeC->config, numofDOFs * sizeof(double));
						TreeC = TreeC -> parent;
					}
					*planlength = N;
					
					//****************//
					return;
				}
			}
			SWAP(&TreeA, &TreeB);
		}
	}
	printf("Fail!!!!!!!!!!!!\n");	
	
	// for(int k = 0; k<N; k++){
		// for(int i=0; i<numofDOFs; i++){
			// printf("%f\t", (*plan)[k][i]);
		// }printf("\n");
	// }
	
	return;
}

//prhs contains input parameters (3): 
//1st is matrix with all the obstacles
//2nd is a row vector of start angles for the arm 
//3nd is a row vector of goal angles for the arm 
//plhs should contain output parameters (2): 
//1st is a 2D matrix plan when each plan[i][j] is the value of jth angle at the ith step of the plan
//(there are D DoF of the arm (that is, D angles). So, j can take values from 0 to D-1
//2nd is planlength (int)
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[])
     
{ 
    
    /* Check for proper number of arguments */    
    if (nrhs != 3) { 
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidNumInputs",
                "Three input arguments required."); 
    } else if (nlhs != 2) {
	    mexErrMsgIdAndTxt( "MATLAB:planner:maxlhs",
                "One output argument required."); 
    } 
        
    /* get the dimensions of the map and the map matrix itself*/     
    int x_size = (int) mxGetM(MAP_IN);
    int y_size = (int) mxGetN(MAP_IN);
    double* map = mxGetPr(MAP_IN);
    
    /* get the start and goal angles*/     
    int numofDOFs = (int) (MAX(mxGetM(ARMSTART_IN), mxGetN(ARMSTART_IN)));
    if(numofDOFs <= 1){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "it should be at least 2");         
    }
    double* armstart_anglesV_rad = mxGetPr(ARMSTART_IN);
    if (numofDOFs != MAX(mxGetM(ARMGOAL_IN), mxGetN(ARMGOAL_IN))){
        	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "numofDOFs in startangles is different from goalangles");         
    }
    double* armgoal_anglesV_rad = mxGetPr(ARMGOAL_IN);
        
    //call the planner
    double** plan = NULL;
    int planlength = 0;
    
    planner(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength); 
    
    printf("planner returned plan of length=%d\n", planlength); 
    
    /* Create return values */
    if(planlength > 0)
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)planlength, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);        
        //copy the values
        int i,j;
        for(i = 0; i < planlength; i++)
        {
            for (j = 0; j < numofDOFs; j++)
            {
                plan_out[j*planlength + i] = plan[i][j];
            }
        }
    }
    else
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);
        //copy the values
        int j;
        for(j = 0; j < numofDOFs; j++)
        {
                plan_out[j] = armstart_anglesV_rad[j];
        }     
    }
    PLANLENGTH_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)1, mxINT8_CLASS, mxREAL); 
    int* planlength_out = (int*)mxGetPr(PLANLENGTH_OUT);
    *planlength_out = planlength;

    
    return;
    
}





