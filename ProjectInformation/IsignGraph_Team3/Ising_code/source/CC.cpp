#include <iostream>
#include <math.h>
#include <stack>
using namespace std;

#define Lx  16
#define Ly  16
#define N  Lx*Ly
#define Three 3


#define max(a,b) (a>b?a:b)
#define min(a,b) (a>b?b:a)

inline int mod(int x, int n)
{
  return  (n + (x % n))%n;
}


inline int nn(int site,int mu)
{
  int x,y;
  int xp, xm, yp, ym;
  int neighbor;
  
  x = site%Lx; y = site/Lx;
  xp = mod(x+1,Lx); xm =  mod(x-1,Lx); yp  = mod(y+1,Ly); ym =  mod(y-1,Ly);
  
  switch(mu){
  case 0:  neighbor =  xp + Lx*y;break; // positive
  case 1:  neighbor =  x + Lx*yp;  break;
  case 2:  neighbor =  xp + Lx*yp; break;
    
  case 3:  neighbor =  xm + Lx*y; break;  // negative
  case 4:  neighbor =  x + Lx*ym; break;
  case 5:  neighbor =  xm + Lx*ym; break;
    
  default: neighbor = -1;
  }
  return neighbor;
}

int  FindGraph(int **frozen,int  *spin, double  *K);
int  FindClusters(int * label, int **frozen, int* size);
void printArray(int * array);

int main()
{

  int  spin[N];
  double K[3];
  K[0] = 4.0;  K[1] = 4.0;  K[2] = 0.0;
  // initialize spins 
  srand(137);
  // Sample Spin Configuration
  for(int s1 = 0; s1 < N ; s1++) rand()%2 == 0 ? spin[s1] = 1 : spin[s1] = -1;
  
  // Define perculation graph
  
  int** frozen = new int*[N];  // row count 
  for(int s = 0; s < N; ++s)    
    frozen[s] = new int[6];   // column count
  
 FindGraph(frozen, spin,K);

#if 0 
  int temp[N];
  for(int mu = 0; mu < 6; mu++)
    {
      cout << "Printf rozen["<< mu<<"] ";
      for(int i = 0; i < N; i++) temp[i] = frozen[i][mu];
      printArray(temp);
    }
#endif
  
 int label[N];
  int size[N];
  int ClusterTotal;
  ClusterTotal = FindClusters(label, frozen, size);

  cout << "spin array  with N = " << N<< endl;
  printArray(spin);
   cout << "Cluster Labels " << endl;
  printArray(label);
  int totSize = 0;
  cout <<"  Total Number of Clusters " <<  ClusterTotal << endl;

  #if 1
  for(int i = 0; i < ClusterTotal ; i++)
    { totSize += size[i];
    cout << "i =  " << i << "   size[i] = "  << size[i] << endl;
    }
  cout<< " Total of Size = " << totSize  << endl;
  #endif
  
  return 0;
}


 // freezing  bonds
int FindGraph(int **frozen, int *spin, double *K)
{
 for(int s1 = 0; s1 < N ; s1++)
   {
     for(int mu = 0; mu <3; mu++)
       {
	 if( spin[nn(s1,mu)]  == spin[s1] &&  (double)rand()/(double)RAND_MAX >  exp(-2*K[mu]) )
	   {
	     frozen[s1][mu] = 1; frozen[nn(s1,mu)][mu+3] = 1;
	   }
	 else
	   {
	     frozen[s1][mu] = 0; frozen[nn(s1,mu)][mu+3] = 0;
	   }
       }
   }
 return 0;
}

// Add int * spin to this rountine if you want it to do  SW update.
int  FindClusters(int* label, int ** frozen,int  * size)
{
 
  for(int site = 0; site < N; site++)
    {
      label[site] = 0; size[site] = 0;
    }  // zero label means it has not been found.
  
  // create the stack
  std::stack<int> stack;
  
  //intialize
  int  seed = 0;
  int  cluster_number = 0;
  
  while(seed < N)
    {
      cluster_number += 1;
      label[seed] = cluster_number;
      stack.push(seed);
      
      // Grow cluster
      while(!stack.empty())
	{
	  int  s_cluster = stack.top();
	  stack.pop();
	  size[cluster_number -1] += 1;
	  for (int mu = 0; mu < 6 ; mu++)
	    {
	     int  site_mu = nn(s_cluster, mu);
	   
	     if( (frozen[s_cluster][mu] == 0) || ( label[site_mu] > 0) )   continue;  // is already found

	    label[site_mu] = cluster_number;
	    stack.push(site_mu);
	  } 
	}  
      while(label[seed] != 0) seed +=1;  // Get next cluster seed
      
    }
  
  return cluster_number; // signed cluster
}


void printArray(int * lattice)
{
	cout<<"\n--------------------------------------------";
	for(int y = 0; y<Ly; y++)
	  {
	     cout << endl;
	    for(int x= 0 ; x<Lx; x++)
   	      printf(" %4d", lattice[x + y* Lx]);    
	  }
	cout<<"\n-------------------------------------------- \n";
}
