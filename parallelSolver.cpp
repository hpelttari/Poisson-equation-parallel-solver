#include <iostream>
#include <stdlib.h>
#include <vector>
#include "mpi.h"

using namespace std;
int const size=6;

//function prototypes
void createGrid(int, double grid[][size]);
void printGrid(double grid[][size]);
void askBoundaryConditions(double grid[][size]);
void updateGridValues(double grid[][size], int, int);
void calculateRedDots(double grid[][size],int, int);
void calculateBlackDots(double grid[][size], int, int);
void updateRedDots(double newGrid[][size],int, int);
void updateBlackDots(double newGrid[][size],int, int);
void sendDots(double newGrid[][size], int);
void receiveDots(double newGrid[][size], double temp[][size], int);
void updateRedDotsOnRootProcess(double newGrid[][size], double temp[][size], int, int);
void updateBlackDotsOnRootProcess(double newGrid[][size], double temp[][size], int, int);


int main(int argc,char ** argv){
  int id,ntasks;

  double grid[size][size];
  
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&id);
  MPI_Comm_size(MPI_COMM_WORLD,&ntasks);

  
  if(id==0){
    createGrid(size,grid);

    cout<<"initial print"<<endl;
    cout<<"\nID="<<id<<endl;
    printGrid(grid);
    cout<<"\n";
    askBoundaryConditions(grid);
  }
  MPI_Bcast(grid,size*size,MPI_DOUBLE,0,MPI_COMM_WORLD);
  /*   cout<<"\nID="<<id<<endl;
       printGrid(grid);
       cout<<"\n";*
  */
  //cout<<"initial grid"<<endl;
  //printGrid(grid);
  for(int i=0;i<100000;i++){
    updateGridValues(grid,id,ntasks);
  }
  /*
    cout<<"final grid"<<endl;
    cout<<"\nID="<<id<<endl;
    printGrid(grid);
    cout<<"\n";
  */
  if(id==0){
    printGrid(grid);
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

void printGrid(double grid[][size]){
  int N=size;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){

      cout<<grid[i][j]<<"         ";
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


void updateGridValues(double newGrid[][size], int id, int ntasks){

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
  calculateRedDots(newGrid, id, ntasks);
  updateRedDots(newGrid,id, ntasks);

  

  MPI_Bcast(newGrid,size*size,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  calculateBlackDots(newGrid, id, ntasks);
  updateBlackDots(newGrid,id,ntasks);
    
  MPI_Bcast(newGrid,size*size,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void calculateRedDots(double newGrid[][size], int id, int ntasks){

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
  //cout<<"red:"<<"id="<<id<<", n="<<n<<", N="<<N<<endl;
  double gamma=1;
  for(int i=n;i<N;i++){
    //cout<<id<<endl;
    //cout<<"id="<<id<<", i="<<i<<", n="<<n<<", N="<<N<<endl;
    //go through every row (every i) every other column (j jumps by two and starts from 1 or 2, depending on the row
    if(i%2==0){
      //cout<<id<<endl;
      for(int j=2;j<size-1;j+=2){
	//cout<<"ida="<<id<<", old:"<<newGrid[i][j]<<endl;
	newGrid[i][j]=(1-gamma)*oldGrid[i][j]+(gamma/4)*(oldGrid[i+1][j]+newGrid[i-1][j]+oldGrid[i][j+1]+newGrid[i][j-1]);
	//cout<<"ida="<<id<<", new:"<<newGrid[i][j]<<endl;
      }
    } else{
      for(int j=1;j<size-1;j+=2){
	//cout<<"idb="<<id<<", old:"<<newGrid[i][j]<<endl;
	newGrid[i][j]=(1-gamma)*oldGrid[i][j]+(gamma/4)*(oldGrid[i+1][j]+newGrid[i-1][j]+oldGrid[i][j+1]+newGrid[i][j-1]);
	//cout<<"idb="<<id<<", new:"<<newGrid[i][j]<<endl;
      }
    }
  }

}

void calculateBlackDots(double newGrid[][size], int id, int ntasks){

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
    //cout<<"black:"<<"id="<<id<<", N="<<N<<", n="<<n<<endl;
  double gamma=1;
  for(int i=n;i<N;i++){
    //cout<<id<<endl;
    //go through every row (every i) every other column (j jumps by two and starts from 1 or 2, depending on the row
    if(i%2==0){
      for(int j=1;j<size-1;j+=2){
	newGrid[i][j]=(1-gamma)*oldGrid[i][j]+(gamma/4)*(oldGrid[i+1][j]+newGrid[i-1][j]+oldGrid[i][j+1]+newGrid[i][j-1]);
      }
    }else{
      for(int j=2;j<size-1;j+=2){
	newGrid[i][j]=(1-gamma)*oldGrid[i][j]+(gamma/4)*(oldGrid[i+1][j]+newGrid[i-1][j]+oldGrid[i][j+1]+newGrid[i][j-1]);
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
  /* int N=(senderID+1)*size/ntasks;
     int n=1+(senderID*size/ntasks);*/


   int N=(senderID+1)*size/ntasks;
  int n=(senderID*size/ntasks);
   if(senderID==ntasks-1){
    N-=1;
  }
  
  for(int i=n;i<N;i++){
    //	cout<<i<<endl;
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
  /* int N=(senderID+1)*size/ntasks;
     int n=1+(senderID*size/ntasks);*/

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


