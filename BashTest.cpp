#include <iostream>
#include <fstream>
#include <omp.h> 
#include <chrono>
#include <iostream>
#include <vector>
#include <set>
using namespace std;
using namespace std::chrono;
#define MAX(a, b) ((a) > (b) ? (a) : (b) )
#define MIN(a, b) ((a) < (b) ? (a) : (b) )
#define White 0
#define Gray 1
#define Black 2
#define INF (1<<31)-1 

int main(int argc, char *argv[])  
{
    const char* filename = (argv[1]);
    int threads = atoi(argv[2]);
	cout<<"\n\n-----------------------------------------------------------------------------------------------\n";
	std::cout << filename<<" threads:"<<threads<<endl;

	return 0;
}