#include <iostream>
#include <fstream>
#include <omp.h> 
#include <chrono>
#include <iostream>
#include <vector>
#include <stack>
#include <set>
#include <bitset>
#include"MyStack.h"
using namespace std;
using namespace std::chrono;
#define MAX(a, b) ((a) > (b) ? (a) : (b) )
#define MIN(a, b) ((a) < (b) ? (a) : (b) )
#define White 0
#define Gray 1
#define Black 2
#define INF (1<<31)-1 
#define bitMapSize 10000000
#define StackSize 1000/// even in large social  network(Twitter), the dinamiter of graph is only 9, 100 should be enouth
//defind node 
typedef struct Vertex {
	int value;
	int source;   // same function as stack defuallt is -1
	int color;
	vector<int> neighbors;    //Vertex's neighbor
	Vertex() : value(0) {}
	//Vertex() : source(-1) {}
	Vertex(int nid) : value(nid) {}
} Vertex;

//Define Graph's adjacency list representation
typedef struct Graph {
	vector<Vertex> vertexs;   //vertex id, same as vertex.value
	int nVertexs;		      //total of node
	bool isDAG;               //wether is direct graph
	int NumDiscover;
	//vector<int> color;        // color
	Graph(int n, bool isDAG) : nVertexs(n), isDAG(isDAG) { vertexs.resize(n); }
	//addEdge
	bool addEdge(int id1, int id2) {
		if (!(MAX(id1, id2) < vertexs.size())) return false;
		if (isDAG) {
			vertexs[id1].neighbors.push_back(id2);
			vertexs[id1].value = id1;
			vertexs[id2].value = id2;
		}
		else {
			vertexs[id1].neighbors.push_back(id2);
			vertexs[id2].neighbors.push_back(id1);
			vertexs[id1].value = id1;
			vertexs[id2].value = id2;
		}
		return true;
	}
	//DFS for demmo 
	vector<int> DFS(int start) {
		set<int> visited;
		vector<int> g, rst;
		g.push_back(start);
		//cout << "push " << start << " ";
		visited.insert(start);
		rst.push_back(start);
		bool found;
		while (g.size() > 0) {
			int id = g[g.size() - 1];
			found = false;
			for (int i = 0; i < vertexs[id].neighbors.size(); i++) {
				int id1 = vertexs[id].neighbors[i];
				if (visited.count(id1) == 0) {
					g.push_back(id1);
					rst.push_back(id1);
					visited.insert(id1);
					//cout << "push " << id1 << " ";
					found = true;
					break;
				}
			}
			if (!found) {
				int id2 = g[g.size() - 1];
				rst.push_back(-1 * id2);
				//cout << "pop " << id2 << " ";
				g.pop_back();
			}
		}
		//cout << endl;
		return rst;
	}
} Graph;
void print_graph(Graph g) {
	for (int i = 0; i < MIN(15, g.nVertexs); i++)
	{
		Vertex curr_vertex = g.vertexs[i];
		int numOfNeighbor = curr_vertex.neighbors.size();

		cout << "node: " << curr_vertex.value << " color: " << curr_vertex.color << " has " << numOfNeighbor << " neighbors ->";
		for (int j = 0; j < numOfNeighbor; j++)
		{
			cout << g.vertexs[i].neighbors[j] << " ";
		}
		cout << "\n";
	}
}
void initialize_vertex(Graph& g) {
	cout << "initialize_vertex...";
	g.NumDiscover = 0;
	for (int i = 0; i < g.nVertexs; i++)
	{

		g.vertexs[i].source = -1;
		g.vertexs[i].color = 0;
		Vertex curr_vertex = g.vertexs[i];
		//cout << "initialize_vertex: " << curr_vertex.value << " color: " << curr_vertex.color << " has " << curr_vertex.neighbors.size() << "\n";
	}
	cout << "finish\n";
}

Graph load_graph(const char* filename) {

	std::ifstream data(filename, std::ios_base::in);
	// email-Eu-core.txt  toy.dat test-1024.txt Wiki-Vote.txt test-5.5kEdges.txt
	if (!data.is_open()) {
		cout << "err: ! file.is_open()";

	}
	int nodeA, nodeB; int num_vertexs = 0;
	//cout << "readdata: "; 
	while (data >> nodeA >> nodeB) {
		num_vertexs = MAX(num_vertexs, MAX(nodeA, nodeB));
	}
	//cout << num_vertexs+1 << " MAX node ID\n ";
	num_vertexs++;//start at 0

	Graph G(num_vertexs, true);
	std::ifstream data2(filename, std::ios_base::in);
	while (data2 >> nodeA >> nodeB) {

		G.addEdge(nodeA, nodeB);
		//cout << nodeA <<"->"<< nodeB << endl;
		//cout << G.vertexs[nodeA].neighbors[0] << " ";
	}
	//cout << "readdata:  finish\n ";
	//print_graph(G);
	return G;
}

void set_color(Graph& g, Vertex v, int color) {
	int id = v.value;
	g.vertexs[id].color = color;
}

void DFS_visit_test(Graph& g, Vertex start) {
	if (start.color == White)
	{
		g.NumDiscover++;
	}
	if (start.color == Black)
	{
		cout << "is Black";
		return;
	}
	//std::vector<Vertex> stack;
	stack<Vertex> stack;
	start.color = 1;
	set_color(g, start, 1);
	stack.push(start);
	while (!stack.empty())
	{
		//Vertex curr = stack.back();
		Vertex curr = stack.top();
		vector<int> neighbors = curr.neighbors;
		//stack.pop_back();
		stack.pop();
		//cout << "thread: " << omp_get_thread_num() << "--";
		cout << "\nat: " << curr.value;           //////////////////////////////////
		for (int i = 0; i < neighbors.size(); i++)
		{
			Vertex u = g.vertexs[neighbors[i]];
			if (i == neighbors.size() - 1 && u.color == Gray)                ///////////////////////////////////////// never  called on toy graph
			{
				curr.color = Black;
				set_color(g, curr, Black);
				cout << " Black: " << curr.value;
			}
			if (u.color == White)
			{
				cout << " discover: " << u.value;                  /////////////////////////////////////////
				u.color = Gray;
				set_color(g, u, 1);
				stack.push(u);
				g.NumDiscover++;
				continue; // if continue search all neighbors, became BFS 
			}
		}
		cout << "\n";                                  /////////////////////////////////////////
	}
}

void DFS_visit_continue(Graph& g, Vertex start) {  //7/10/2020 not use new vector
	if (start.color == White)
	{
		g.NumDiscover++;
	}
	if (start.color == Black)
	{
		//cout << "is Black";
		return;
	}
	//std::vector<Vertex> stack;
	stack<int> stack;
	//start.color = 1;
	//set_color(g, start, 1);
	g.vertexs[start.value].color = Gray;
	stack.push(start.value);
	while (!stack.empty())
	{
		//Vertex curr = stack.back();
		int currValue = stack.top();
		int currColor = g.vertexs[currValue].color;
		stack.pop();

		if (g.vertexs[currValue].neighbors.size()==0)
		{
			g.vertexs[currValue].color = Black;

			continue;
		}
		if (g.vertexs[currValue].neighbors.size() == 1)
		{
			//int oneneighborValue = g.vertexs[currValue].neighbors[0];
			g.vertexs[currValue].color = Black;
			stack.push(g.vertexs[currValue].neighbors[0]);
			continue;
		}
		//if (g.vertexs[currValue].neighbors.size() == 2)
		//{


		//	continue;
		//}
		// if neighbors lager than 3 cause new vector<int> is expansive
		//vector<int> neighbors = g.vertexs[currValue].neighbors;
		
		
		//cout << "thread: " << omp_get_thread_num() << "--";
		//cout << "\nat: " << curr.value;           //////////////////////////////////
		for (int i = 0; i < g.vertexs[currValue].neighbors.size(); i++)
		{
			//Vertex u = g.vertexs[neighbors[i]];
			int Ucolor = g.vertexs[g.vertexs[currValue].neighbors[i]].color;
			int Uvalue = g.vertexs[g.vertexs[currValue].neighbors[i]].value;
			//if (i == neighbors.size() - 1 && u.color == Gray)                ///////////////////////////////////////// never  called on toy graph
			if (i == g.vertexs[currValue].neighbors.size() - 1 && Ucolor == Gray)
			{
				currColor = Black;
				//set_color(g, curr, Black);
				g.vertexs[currValue].color = Black;
				//cout << " Black: " << curr.value;
			}
			if (Ucolor == White)
			{
				//cout << " discover: " << u.value;                  /////////////////////////////////////////
				Ucolor = Gray;
				//set_color(g, u, 1);
				g.vertexs[Uvalue].color = Gray;
				stack.push(Uvalue);
				g.NumDiscover++;
				continue; // if continue search all neighbors, became  DFS break
			}
		}
		//cout << "\n";                                  /////////////////////////////////////////
	}
}
void DFS_visit(Graph& g, Vertex start) {
	if (start.color == White)
	{
		g.NumDiscover++;
	}
	if (start.color == Black)
	{
		return;
	}
	std::vector<Vertex> stack;
	start.color = 1;
	set_color(g, start, 1);
	stack.push_back(start);
	while (!stack.empty())
	{
		Vertex curr = stack.back();
		vector<int> neighbors = curr.neighbors;
		stack.pop_back();
		for (int i = 0; i < neighbors.size(); i++)
		{
			Vertex u = g.vertexs[neighbors[i]];
			if (i == neighbors.size())                ///////////////////////////////////////// never  called on toy graph
			{
				curr.color = Black;
				set_color(g, curr, Black);
			}
			if (u.color == White)
			{
				u.color = Gray;
				set_color(g, u, 1);
				stack.push_back(u);
				g.NumDiscover++;
				break; // if continue search all neighbors, became BFS  break
			}
		}                             /////////////////////////////////////////
	}

}

void DFS_visit_continue_bitmap(Graph& g, Vertex start, bitset<bitMapSize>& color_bit) {  //7/13/2020 not use new vector
	int startId = start.value;
	if (!color_bit.test(start.value)) //// if node i's color not 1 : !color_bit.test(i)
	{
		g.NumDiscover++;
	}
	else
	{
		return;
	} 
	register stack<int> stack;

	//g.vertexs[start.value].color = Gray;
	color_bit.set(start.value,1);

	stack.push(start.value);
	while (!stack.empty())
	{
		int currValue = stack.top();
		//int currColor = g.vertexs[currValue].color;
		int currColor = color_bit.test(currValue) ? Gray : White;/// if ture, currcolor ==1(gray)
		stack.pop();

        
		if (g.vertexs[currValue].neighbors.size()==0)
		{
			// g.vertexs[currValue].color = Black;
            color_bit.set(currValue,1);
			continue;
		}
		if (g.vertexs[currValue].neighbors.size() == 1)
		{
			//int oneneighborValue = g.vertexs[currValue].neighbors[0];
			// g.vertexs[currValue].color = Black;
            color_bit.set(currValue,1);
			stack.push(g.vertexs[currValue].neighbors[0]);
			continue;
		}


        

		for (int i = 0; i < g.vertexs[currValue].neighbors.size(); i++)
		{
			//Vertex u = g.vertexs[neighbors[i]];
			//int Ucolor = g.vertexs[g.vertexs[currValue].neighbors[i]].color;
			int Ucolor = color_bit.test(g.vertexs[currValue].neighbors[i]) ? Gray : White;
			int Uvalue = g.vertexs[g.vertexs[currValue].neighbors[i]].value;
	
			//if (i == g.vertexs[currValue].neighbors.size() - 1 && Ucolor == Gray)////////////// no black now
			//{
			//	currColor = Black;
			//	//set_color(g, curr, Black);
			//	g.vertexs[currValue].color = Black;
			//	//cout << " Black: " << curr.value;
			//}

			if (Ucolor == White)
			{
				//cout << " discover: " << u.value;                  /////////////////////////////////////////

				//g.vertexs[Uvalue].color = Gray;
				color_bit.set(Uvalue, Gray);

				stack.push(Uvalue);
				g.NumDiscover++;
				continue; // if continue search all neighbors, became  DFS break
			}
		}
		//cout << "\n";                                  /////////////////////////////////////////
	}
}


void DFS_visit_continue_removeSTLstack(Graph& g, Vertex start) {  //   //7/14/2020 use static stack instead of STL Stack is my, stack is STL
	Stack<int> stack(1000);
	if (start.color == White)
	{
		g.NumDiscover++;
	}
	if (start.color == Black)
	{
		return;
	}

	g.vertexs[start.value].color = Gray;
	stack.push(start.value);
	while (!stack.isEmpty())
	{
		int currValue = stack.getTop();
		int currColor = g.vertexs[currValue].color;
		int x = 0;
		stack.pop(&x);


		for (int i = 0; i < g.vertexs[currValue].neighbors.size(); i++)
		{
			//Vertex u = g.vertexs[neighbors[i]];
			int Ucolor = g.vertexs[g.vertexs[currValue].neighbors[i]].color;
			int Uvalue = g.vertexs[g.vertexs[currValue].neighbors[i]].value;
			//if (i == neighbors.size() - 1 && u.color == Gray)                ///////////////////////////////////////// never  called on toy graph
			if (i == g.vertexs[currValue].neighbors.size() - 1 && Ucolor == Gray)
			{
				currColor = Black;
				//set_color(g, curr, Black);
				g.vertexs[currValue].color = Black;
				//cout << " Black: " << curr.value;
			}
			if (Ucolor == White)
			{
				//cout << " discover: " << u.value;                  /////////////////////////////////////////
				Ucolor = Gray;
				//set_color(g, u, 1);
				g.vertexs[Uvalue].color = Gray;
				stack.push(Uvalue);
				g.NumDiscover++;
				continue; // if continue search all neighbors, became  DFS break
			}
		}
		//cout << "\n";                                  /////////////////////////////////////////
	}
}

void DFS_visit_removeGraph(int** Adjacency, int start, int& NumDiscover) {  //   //7/15/2020 removeGraph from DFS_visit_continue_removeSTLstack
	Stack<int> stack(1000);
	
	if (Adjacency[start][0] == White)
	{
		NumDiscover++;
	}
	if (Adjacency[start][0] == Black)
	{
		return;
	}

	Adjacency[start][0] = Gray;
	stack.push(start);
	while (!stack.isEmpty())
	{
		int currValue = 0;
		stack.pop(&currValue);
		int currColor = Adjacency[start][0];
		
		                ///////////////////////////////////////// ///////////////////////////////////////// TODO pop and push has problem
		//cout << currValue<< "poped, stack.size() " << stack.size() << endl;
		//cout << " at: " << currValue;
		if (Adjacency[currValue][1]==0)
		{
			continue;
		}
		for (int i = 2; i < Adjacency[currValue][1]+2; i++)
		{
			//Vertex u = g.vertexs[neighbors[i]];
			int Uvalue = Adjacency[currValue][i];
			int Ucolor = Adjacency[Uvalue][0];
			//int Uvalue = g.vertexs[g.vertexs[currValue].neighbors[i]].value;
			//if (i == neighbors.size() - 1 && u.color == Gray)                ///////////////////////////////////////// never  called on toy graph
			if (i == Adjacency[currValue][1] - 1 && Ucolor == Gray)
			{
				//g.vertexs[currValue].color = Black;
				Adjacency[currValue][0] = Black;
				//cout << " Black: " << currValue;              /////////////////////////////////////////
			}
			if (Ucolor == White)
			{
				//cout << " discover: " << Uvalue;                  /////////////////////////////////////////
				//Ucolor = Gray;
				//set_color(g, u, 1);
				Adjacency[Uvalue][0]  = Gray;
				stack.push(Uvalue);
				NumDiscover++;
				continue; // if continue search all neighbors, became  DFS break
			}
		}
		//cout << "\n";                                  /////////////////////////////////////////
	}
}

void DFS_test(Graph& g, int start) {
	for (int i = start; i < g.nVertexs; i++)
	{
		if (g.NumDiscover >= g.nVertexs)
		{
			cout << "g.NumDiscover: " << g.NumDiscover << endl;
			return;
		}
		Vertex curr = g.vertexs[i];
		//cout << curr.color;
		if (curr.color == White)
		{
			DFS_visit_test(g, curr);
		}
	}
	cout << "g.NumDiscover: " << g.NumDiscover << endl;
}
void DFS(Graph& g, int start) {
	for (int i = start; i < g.nVertexs; i++)
	{

		//Vertex curr = g.vertexs[i];
		//cout << curr.color;
		if (g.vertexs[i].color == White)
		{
			DFS_visit(g, g.vertexs[i]);
		}

		if (g.NumDiscover >= g.nVertexs)
		{
			//cout << "g.NumDiscover: " << g.NumDiscover  << endl;
			return;
		}
	}
	//cout << "g.NumDiscover: " << g.NumDiscover << endl;
}
void DFS_optimize2(Graph& g, int start) { // use DFS_visit_continue to change DFS_visit.  sequencial ordered
	for (int i = start; i < g.nVertexs; i++)
	{

		//Vertex curr = g.vertexs[i];
		//cout << curr.color;
		if (g.vertexs[i].color == White)
		{
			DFS_visit_continue(g, g.vertexs[i]);
		}

		if (g.NumDiscover >= g.nVertexs)
		{
			//cout << "g.NumDiscover: " << g.NumDiscover  << endl;
			return;
		}
	}
	//cout << "g.NumDiscover: " << g.NumDiscover << endl;
}
void DFS_optimize3_bitmap(Graph& g, int start, bitset<bitMapSize>& color_bit) { // use bit map to store color info 0 is flase(unvisited), 1 is true, visited
	for (int i = start; i < g.nVertexs; i++)
	{
		if (!color_bit.test(i))// if node i's color not 1 : !color_bit.test(i)
		{
			DFS_visit_continue_bitmap(g, g.vertexs[i], color_bit);
		}

		if (color_bit.count() >= g.nVertexs)
		{
			return;
		}
	}
}
void DFS_optimize4_removeSTLstack(Graph& g, int start) { //7/14/2020 use static stack instead of STL 
	for (int i = start; i < g.nVertexs; i++)
	{
		if (g.vertexs[i].color == White)
		{
			DFS_visit_continue_removeSTLstack(g, g.vertexs[i]);
		}
		if (g.NumDiscover >= g.nVertexs)
		{return;
		}
	}//cout << "g.NumDiscover: " << g.NumDiscover << endl;
}
void DFS_optimize5_removeGraphStruct(int** Adjacency, int start,int numofVertex) { //7/15/2020 removeGraphStruct, other same as optimize4
	int NumDiscover = 0;//////////use int* in mutil thread or use 
	for (int i = start; i < numofVertex; i++)
	{
		if (Adjacency[i][0] == White)
		{
			DFS_visit_removeGraph(Adjacency, i, NumDiscover);
		}
		if (NumDiscover >= numofVertex)
		{
			return;
		}
	}//cout << "g.NumDiscover: " << g.NumDiscover << endl;
}


//put current visiting state of graph color info in to a bit map to save memeory
void colorBitmap(Graph g, bitset<bitMapSize> color_bit) {
	cout << "colorBitmap...";
	g.NumDiscover = 0;
	for (int i = 0; i < g.nVertexs; i++)
	{
		color_bit.set(i, g.vertexs[i].color);
		//cout << "initialize_vertex: " << curr_vertex.value << " color: " << curr_vertex.color << " has " << curr_vertex.neighbors.size() << "\n";
	}
	cout << "finish\n";

}


void printlist(int* a, int size) {
	for (int i = 0; i < size; i++)
	{
		cout << a[i] << " ";
	}cout << endl;
}

void CreateAdjacency( int** Adjacency, int numofVertex, Graph g) {
	
	int end = numofVertex;
	cout <<end<< " size CreateAdjacency...\n";
	for (int i = 0; i < end; i++)
	{
		Vertex curr = g.vertexs[i];
		vector<int> neighbors = curr.neighbors;
		int numofneighbor = neighbors.size();

		int size = numofneighbor + 2;
		Adjacency[i] = new int[size];
		Adjacency[i][0] = 0; Adjacency[i][1] = numofneighbor;  //default color is 0
		//cout << "siez of Adjacency[i] : " << sizeof(Adjacency[i]) / sizeof(Adjacency[i][0])<<endl;
		if (numofneighbor>0)
		{
			for (int j = 2; j < size; j++)
			{
				Adjacency[i][j] = neighbors[ j-2];
			}
		}

		//printlist(Adjacency[i], size);
	}
	cout << "finish\n";
}


int main(int argc, char *argv[])  
{
    const char* filename = (argv[1]);
    // int threads = atoi(argv[2]);
	register std::bitset<bitMapSize> color_bit;
	color_bit.reset();
	//std::bitset<1000> foo(std::string("1011"));
	//std::cout << foo.reset() << '\n';
	//std::cout << "foo.size() is " << foo.size()/(1024*1024) << '\n';
	std::cout << "size of bitmap is(MB): " << 1.0*sizeof(color_bit) / (1024 * 1024) << '\n';

	auto start = high_resolution_clock::now();
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);




	int threads = 2;
	// const char* filename = "com-orkut.ungraph.txt";
	cout << " filename " << filename;
	Graph g = load_graph(filename); //
	// email-Eu-core.txt  toy.dat test-1024.txt Wiki-Vote.txt test-5.5kEdges.txt   p2p-Gnutella31.txt
	// toy_DCG.dat  com-orkut.ungraph.txt

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	initialize_vertex(g);
	int numofVertex = g.nVertexs;
	int** Adjacency = new int* [numofVertex];
	CreateAdjacency(Adjacency, numofVertex, g);
	start = high_resolution_clock::now();
	//////////////////////////////////
	DFS_optimize5_removeGraphStruct(Adjacency, 0, numofVertex);
	//////////////////////////////////
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	std::cout << "\n DFS_optimize5_removeGraphStruct time spends (ms): " << duration.count() << endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// initialize_vertex(g); // set all color as 0, sources as -1
	// start = high_resolution_clock::now();
	// //////////////////////////////////////'
	// DFS_optimize3_bitmap(g, 0, color_bit);

	// //////////////////////////////////
	// stop = high_resolution_clock::now();
	// duration = duration_cast<milliseconds>(stop - start);
	// std::cout << "\n DFS_optimize3_bitmap  serial time spends (ms): " << duration.count() << endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	initialize_vertex(g); // set all color as 0, sources as -1
	start = high_resolution_clock::now();
	//////////////////////////////////////'
	DFS_optimize2(g, 0);

	//////////////////////////////////
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	std::cout << "\n nDFS_optimize2  serial time spends (ms): " << duration.count() << endl;

	initialize_vertex(g);
	start = high_resolution_clock::now();
	//////////////////////////////////////'
	DFS(g, 0);
	//////////////////////////////////
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	std::cout << "\n serial time spends (ms): " << duration.count() << endl;

	
	
	
	
	return 0;
}