#include <iostream>
#include <vector>

using namespace std;
int const size=10;

//function prototypes
void createGrid(int, double grid[][size]);
void printGrid(double grid[][size]);
void askBoundaryConditions(double grid[][size]);
void updateGridValues(double grid[][size]);
void calculateRedDots(double grid[][size]);
void calculateBlackDots(double grid[][size]);


int main(){
  double grid[size][size];
  createGrid(size,grid);
  printGrid(grid);
  askBoundaryConditions(grid);
  cout<<"initial grid"<<endl;
  printGrid(grid);
  for(int i=0;i<10000000;i++){
    updateGridValues(grid);
  }
  cout<<"new grid"<<endl;
  printGrid(grid);
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
}


void updateGridValues(double newGrid[][size]){

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
  calculateRedDots(newGrid);
  calculateBlackDots(newGrid);
}

void calculateRedDots(double newGrid[][size]){

  int N=size;
  double oldGrid[size][size];
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      oldGrid[i][j]=newGrid[i][j];
    }
  }

  
  double gamma=1;
  for(int i=1;i<N-1;i++){
    //go through every row (every i) every other column (j jumps by two and starts from 1 or 2, depending on the row
    if(i%2==0){
    for(int j=2;j<N-1;j+=2){
      newGrid[i][j]=(1-gamma)*oldGrid[i][j]+(gamma/4)*(oldGrid[i+1][j]+newGrid[i-1][j]+oldGrid[i][j+1]+newGrid[i][j-1]);
    }
    } else{
    for(int j=1;j<N-1;j+=2){
      newGrid[i][j]=(1-gamma)*oldGrid[i][j]+(gamma/4)*(oldGrid[i+1][j]+newGrid[i-1][j]+oldGrid[i][j+1]+newGrid[i][j-1]);
    }
    }
  }
}

void calculateBlackDots(double newGrid[][size]){

  int N=size;
  double oldGrid[size][size];
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      oldGrid[i][j]=newGrid[i][j];
    }
  }
  
  double gamma=1;
  for(int i=1;i<N-1;i++){
    //go through every row (every i) every other column (j jumps by two and starts from 1 or 2, depending on the row
    if(i%2==0){
    for(int j=1;j<N-1;j+=2){
      newGrid[i][j]=(1-gamma)*oldGrid[i][j]+(gamma/4)*(oldGrid[i+1][j]+newGrid[i-1][j]+oldGrid[i][j+1]+newGrid[i][j-1]);
    }
    }else{
    for(int j=2;j<N-1;j+=2){
      newGrid[i][j]=(1-gamma)*oldGrid[i][j]+(gamma/4)*(oldGrid[i+1][j]+newGrid[i-1][j]+oldGrid[i][j+1]+newGrid[i][j-1]);
    }
      
    }
  }
}
