#include <iostream>
#include <stdlib.h>
#include <vector>
#include "mpi.h"
#include <cmath>
#include <ctime>
#include <iomanip>

using namespace std;
//using std::setw;
int const size=15;

//function prototypes
void createGrid(int, double grid[][size]);
void printGrid(double grid[][size]);
void askBoundaryConditions(double grid[][size]);
void updateGridValues(double grid[][size], int, int, double);
void calculateRedDots(double grid[][size],int, int, double);
void calculateBlackDots(double grid[][size], int, int, double);
void updateRedDots(double newGrid[][size],int, int);
void updateBlackDots(double newGrid[][size],int, int);
void sendDots(double newGrid[][size], int);
void receiveDots(double newGrid[][size], double temp[][size], int);
void updateRedDotsOnRootProcess(double newGrid[][size], double temp[][size], int, int);
void updateBlackDotsOnRootProcess(double newGrid[][size], double temp[][size], int, int);
void copyGrid1ToGrid2(double grid1[][size],double grid2[][size]);
void initializeErrorMatrix(double errorMatrix[][size], int);
double g(int,int);
double calculateErrorMatrixAndReturnMaxError(double currentGrid[][size],double previousGrid[][size],double errorGrid[][size]);


int main(int argc,char ** argv){
  int id,ntasks;

  double grid[size][size];
  double previousGrid[size][size];
  double errorGrid[size][size];
  double maxError;
  double gamma=0.95;
  double convergenceCriterion=0.00001;
  clock_t t1,t2;
  bool printOn=0;


  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&id);
  MPI_Comm_size(MPI_COMM_WORLD,&ntasks);

    if(id==0){
    cout<<"press 1 to print grid values, or 0 to disable printing:"<<endl;
    char val;
    cin>>val;
    if(val=='1'){
      printOn=1;
    }
  }

  //check that size is divisible with the number of processors
  if(size%ntasks!=0){
    if(id==0){
      cout<<"The size of the array must be divisible with the number of processors used"<<endl;
      cout<<"The size of the current array is "<<size<<endl;
      cout<<"Run the program again with different number of processors"<<endl;
    }
    return 1;
  }
  
  if(id==0){
    createGrid(size,grid);
    initializeErrorMatrix(errorGrid,size);
    if(printOn){
      cout<<"initial print"<<endl;
      cout<<"\nID="<<id<<endl;
      printGrid(grid);
    }
    cout<<"\n";
    askBoundaryConditions(grid);
    copyGrid1ToGrid2(grid,previousGrid);
    if(printOn){
      cout<<"Boundary values added:"<<endl;
      printGrid(grid);
      cout<<"Initial error Grid:"<<endl;
      printGrid(errorGrid);
      cout<<endl;
    }
  }
  
  MPI_Bcast(grid,size*size,MPI_DOUBLE,0,MPI_COMM_WORLD);
  t1=clock();
  while(true){
    if(id==0){
    copyGrid1ToGrid2(grid,previousGrid);
    }
    
    updateGridValues(grid,id,ntasks,gamma);

    if(id==0){
    maxError=calculateErrorMatrixAndReturnMaxError(grid,previousGrid,errorGrid);
    }
    MPI_Bcast(&maxError,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(maxError<convergenceCriterion){
      break;
    }
  }
  t2=clock();

  if(id==0){
    if(printOn){
      cout<<"Final solved grid:"<<endl;
    printGrid(grid);
    cout<<endl;
    cout<<"errors:"<<endl;
    printGrid(errorGrid);
    }
    cout<<"CPU time: "<<(double)(t2-t1)/CLOCKS_PER_SEC<<endl;
  }
  MPI_Finalize();
  return 0;
}

void createGrid(int N, double grid[][size]){
  double guess;
  cout<<"Give the initial value for interior points:"<<endl;
  cin>>guess;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      grid[i][j]=guess;
    }
  }
}

void initializeErrorMatrix(double errorMatrix[][size], int size){
  for(int i=0;i<size;i++){
    for(int j=0;j<size;j++){
      errorMatrix[i][j]=0;
    }
  }
}

void copyGrid1ToGrid2(double grid1[][size],double grid2[][size]){
  for(int i=0;i<size;i++){
    for(int j=0;j<size;j++){
      grid2[i][j]=grid1[i][j];
    }
  }
  
}

double calculateErrorMatrixAndReturnMaxError(double currentGrid[][size],double previousGrid[][size], double errorGrid[][size]){
  double maxError=0;
  for(int i=1;i<size-1;i++){
    for(int j=1;j<size-1;j++){
      errorGrid[i][j]=abs((currentGrid[i][j]-previousGrid[i][j])/currentGrid[i][j]);
      if(errorGrid[i][j]>maxError){
	maxError=errorGrid[i][j];
      
      }
    }
  }
    return maxError;
}

void printGrid(double grid[][size]){
  int N=size;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){

      cout<<left<<setw(13)<<grid[i][j]<<setw(1)<<" ";
    }
    cout<<"\n";
  }
}


void askBoundaryConditions(double grid[][size]){
  double inputGrid[size][size];
  char whichCondition;
  int N=size;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      inputGrid[i][j]=grid[i][j];
    }
  }

  cout<<"Type a if you want all the sides have the same boundary condition and b if you want different conditions for each side"<<endl;
  cin>>whichCondition;

  if(whichCondition=='a'){
    double edge;
    //boundary conditions same at each edge:
    cout<<"Give boundary conditions"<<endl;
    cout<<"Give the value of f(x,y) at the edges:"<<endl;
    cin>>edge;
    for(int i=0;i<N;i++){
      grid[0][i]=edge;
      grid[i][N-1]=edge;
      grid[i][0]=edge;
      grid[N-1][i]=edge;
    }
  }
  
  //TODO: boundary conditions not same at every edge:

  if(whichCondition=='b'){
    double leftSide,rightSide,lowerSide,upperSide;

    cout<<"Give boundary conditions:"<<endl;
    cout<<"conditon for the left side"<<endl;
    cin>>leftSide;
    cout<<"lower side"<<endl;
    cin>>lowerSide;
    cout<<"upper side"<<endl;
    cin>>upperSide;
    cout<<"right side"<<endl;
    cin>>rightSide;
    for(int i=0;i<N;i++){
      grid[0][i]=upperSide;
      grid[N-1][i]=lowerSide;
    }

    for(int i=0;i<N;i++){
      grid[i][0]=leftSide;
      grid[i][N-1]=rightSide;
    }

  
  }
  cout<<endl;
}


void updateGridValues(double newGrid[][size], int id, int ntasks, double gamma){

  int N=size;
  double oldGrid[size][size];
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      oldGrid[i][j]=newGrid[i][j];
    }
  }
  /*double gamma=1;
    int N = newGrid.size();
    for(int i=1;i<N-1;i++){
    for(int j=1;j<N-1;j++){
    newGrid[i][j]=(1-gamma)*oldGrid[i][j]+(gamma/4)*(oldGrid[i+1][j]+newGrid[i-1][j]+oldGrid[i][j+1]+newGrid[i][j-1]);
    }
    }*/
  calculateRedDots(newGrid, id, ntasks, gamma);
  updateRedDots(newGrid,id, ntasks);

  

  MPI_Bcast(newGrid,size*size,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  calculateBlackDots(newGrid, id, ntasks, gamma);
  updateBlackDots(newGrid,id,ntasks);
    
  MPI_Bcast(newGrid,size*size,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void calculateRedDots(double newGrid[][size], int id, int ntasks, double gamma){

  int N=size;
  double oldGrid[size][size];
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      oldGrid[i][j]=newGrid[i][j];
    }
  }



  N=(id+1)*size/ntasks;
  int n=(id*size/ntasks);
  if(id==0){
    n+=1;
  } else if(id==ntasks-1){
    N-=1;
  }
  //double gamma=1;
  for(int i=n;i<N;i++){
    //go through every row (every i) every other column (j jumps by two and starts from 1 or 2, depending on the row
    if(i%2==0){
      for(int j=2;j<size-1;j+=2){
	newGrid[i][j]=(1-gamma)*oldGrid[i][j]+(gamma/4)*(oldGrid[i+1][j]+newGrid[i-1][j]+oldGrid[i][j+1]+newGrid[i][j-1])-(gamma/(4.0*size))*g(i,j);
      }
    } else{
      for(int j=1;j<size-1;j+=2){
	newGrid[i][j]=(1-gamma)*oldGrid[i][j]+(gamma/4)*(oldGrid[i+1][j]+newGrid[i-1][j]+oldGrid[i][j+1]+newGrid[i][j-1])-(gamma/(4.0*size))*g(i,j);
      }
    }
  }

}

void calculateBlackDots(double newGrid[][size], int id, int ntasks, double gamma){

  int N=size;
  double oldGrid[size][size];
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      oldGrid[i][j]=newGrid[i][j];
    }
  }

  N=(id+1)*size/ntasks;
  int n=(id*size/ntasks);
    if(id==0){
    n+=1;
  } else if(id==ntasks-1){
    N-=1;
  }
    //double gamma=1;
  for(int i=n;i<N;i++){
    //go through every row (every i) every other column (j jumps by two and starts from 1 or 2, depending on the row
    if(i%2==0){
      for(int j=1;j<size-1;j+=2){
	newGrid[i][j]=(1-gamma)*oldGrid[i][j]+(gamma/4)*(oldGrid[i+1][j]+newGrid[i-1][j]+oldGrid[i][j+1]+newGrid[i][j-1])-(gamma/(4.0*size))*g(i,j);
      }
    }else{
      for(int j=2;j<size-1;j+=2){
	newGrid[i][j]=(1-gamma)*oldGrid[i][j]+(gamma/4)*(oldGrid[i+1][j]+newGrid[i-1][j]+oldGrid[i][j+1]+newGrid[i][j-1])-(gamma/(4.0*size))*g(i,j);
      }
      
    }
  }
}

void updateRedDots(double newGrid[][size],int id, int ntasks){
  if(id!=0){
    sendDots(newGrid,id);
  } else{
    double temp[size][size];
    for(int senderID=1;senderID<ntasks;senderID++){
      receiveDots(newGrid, temp, senderID);
      updateRedDotsOnRootProcess(newGrid,temp,senderID,ntasks);
    }
  }
  
}

void updateBlackDots(double newGrid[][size],int id, int ntasks){
  if(id!=0){
    sendDots(newGrid,id);
  } else{
    double temp[size][size];
    for(int senderID=1;senderID<ntasks;senderID++){
      receiveDots(newGrid, temp, senderID);
      updateBlackDotsOnRootProcess(newGrid,temp,senderID,ntasks);
    }
  }
  
}



void sendDots(double newGrid[][size],int id){
  MPI_Send(newGrid,size*size,MPI_DOUBLE,0,id,MPI_COMM_WORLD);
}

void receiveDots(double newGrid[][size], double temp[][size],int senderID){
  MPI_Recv(temp,size*size,MPI_DOUBLE,senderID,senderID,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

void updateRedDotsOnRootProcess(double newGrid[][size], double temp[][size], int senderID, int ntasks){

   int N=(senderID+1)*size/ntasks;
  int n=(senderID*size/ntasks);
   if(senderID==ntasks-1){
    N-=1;
  }
  
  for(int i=n;i<N;i++){
    //go through every row (every i) every other column (j jumps by two and starts from 1 or 2, depending on the row
    if(i%2==0){
      for(int j=2;j<size-1;j+=2){
	newGrid[i][j]=temp[i][j];
      }
    } else{
      for(int j=1;j<size-1;j+=2){
	newGrid[i][j]=temp[i][j];
      }
    }
  }
}

void updateBlackDotsOnRootProcess(double newGrid[][size], double temp[][size], int senderID, int ntasks){

     int N=(senderID+1)*size/ntasks;
  int n=(senderID*size/ntasks);
   if(senderID==ntasks-1){
    N-=1;
  }

  
  for(int i=n;i<N;i++){
    //go through every row (every i) every other column (j jumps by two and starts from 1 or 2, depending on the row
    if(i%2==0){
      for(int j=1;j<size-1;j+=2){
	newGrid[i][j]=temp[i][j];
      }
    } else{
      for(int j=2;j<size-1;j+=2){
	newGrid[i][j]=temp[i][j];
      }
    }
  }
      
}

double g(int i,int j){
  double x=(double)i/size;
  double y=(double)j/size;
  return x*y;
}


