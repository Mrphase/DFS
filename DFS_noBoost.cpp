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
	for (int i = 0; i < MIN(15,g.nVertexs); i++)
	{
		Vertex curr_vertex = g.vertexs[i];
		int numOfNeighbor = curr_vertex.neighbors.size();
		
		cout<< "node: "<< curr_vertex.value<<" color: "<<curr_vertex.color<<" has "<< numOfNeighbor << " neighbors ->";
		for (int j = 0; j < numOfNeighbor; j++)
		{
			cout << g.vertexs[i].neighbors[j]<<" ";
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
Graph load_graph_cout(const char* filename) {

	std::ifstream data(filename, std::ios_base::in);
	// email-Eu-core.txt  toy.dat test-1024.txt Wiki-Vote.txt test-5.5kEdges.txt
	if (!data.is_open()) {
		cout << "err: ! file.is_open()";

	}
	int nodeA, nodeB; int num_vertexs = 0;
	cout << "readdata: ";
	while (data >> nodeA >> nodeB) {
		num_vertexs = MAX(num_vertexs, MAX(nodeA, nodeB));
	}cout << num_vertexs + 1 << " MAX node ID\n "; num_vertexs++;//start at 0

	Graph G(num_vertexs, true);
	std::ifstream data2(filename, std::ios_base::in);
	while (data2 >> nodeA >> nodeB) {

		G.addEdge(nodeA, nodeB);
		//cout << nodeA <<"->"<< nodeB << endl;
		//cout << G.vertexs[nodeA].neighbors[0] << " ";
	}cout << "readdata:  finish\n ";
	//print_graph(G);
	return G;
}
Graph load_graph( const char* filename) {

	std::ifstream data(filename, std::ios_base::in);
	// email-Eu-core.txt  toy.dat test-1024.txt Wiki-Vote.txt test-5.5kEdges.txt
		if (!data.is_open()) {
		cout << "err: ! file.is_open()";
		
	}
		int nodeA, nodeB; int num_vertexs = 0;
	//cout << "readdata: "; 
	while (data >> nodeA >> nodeB) {
		num_vertexs=MAX(num_vertexs,MAX(nodeA,nodeB));
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
	if (start.color==White)
	{
		g.NumDiscover++;
	}
	if (start.color==Black)
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
		//cout << "thread: " << omp_get_thread_num() << "--";
		cout << "\nat: " << curr.value;           //////////////////////////////////
		for (int i = 0; i < neighbors.size(); i++)
		{
			Vertex u = g.vertexs[neighbors[i]];
			if (i == neighbors.size())                ///////////////////////////////////////// never  called on toy graph
			{
				curr.color = Black;
				set_color(g, curr, Black);
				cout << " finish: " << curr.value;
			}
			if (u.color== White)
			{


				cout << " discover: " << u.value;                  /////////////////////////////////////////
				u.color = Gray;
				set_color(g, u, 1);
				stack.push_back(u);
				g.NumDiscover++;
				break; // if continue search all neighbors, became BFS 
			}


		}
		cout << "\n";                                  /////////////////////////////////////////
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
		//cout << "thread: " << omp_get_thread_num() << "--";
		//cout << "\nat: " << curr.value;           //////////////////////////////////
		for (int i = 0; i < neighbors.size(); i++)
		{
			Vertex u = g.vertexs[neighbors[i]];
			if (i == neighbors.size())                ///////////////////////////////////////// never  called on toy graph
			{
				curr.color = Black;
				set_color(g, curr, Black);
				//cout << " finish: " << curr.value;
			}
			if (u.color == White)
			{


				//cout << " discover: " << u.value;                  /////////////////////////////////////////
				u.color = Gray;
				set_color(g, u, 1);
				stack.push_back(u);
				g.NumDiscover++;
				break; // if continue search all neighbors, became BFS 
			}


		}
		//cout << "\n";                                  /////////////////////////////////////////
	}

}



void DFS_visit_workload(Graph& g, Vertex start, int* &workload) {
	if (start.color == White)
	{
		g.NumDiscover++;
		workload[omp_get_thread_num()]++;   //////// record workload
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
				workload[omp_get_thread_num()]++;   //////// record workload

				break; // if continue search all neighbors, became BFS 
			}
		}                             /////////////////////////////////////////
	}

}

void DFS_visit_unOrder(Graph& g, Vertex start,int  threadId, int threads) {
	if (start.color == White)
	{
		g.NumDiscover++;
	}
	if (start.color == Black)
	{
		return;
	}
	/// else color is gray
	std::vector<Vertex> stack;
	start.color = 1;
	set_color(g, start, 1);
	stack.push_back(start);
	while (!stack.empty())
	{
		Vertex curr = stack.back();
		vector<int> neighbors = curr.neighbors;
		stack.pop_back();
		//cout << "thread: " << omp_get_thread_num() << "--";
		//cout << "\nat: " << curr.value;           //////////////////////////////////
		for (int i = threadId; i < neighbors.size(); i+= threads)        //////////////////////////////////    for 6 threads
		{
			Vertex u = g.vertexs[neighbors[i]];
			if (i == neighbors.size())                ///////////////////////////////////////// never  called on toy graph
			{
				curr.color = Black;
				set_color(g, curr, Black);
				cout << " finish: " << curr.value;
			}
			if (u.color == White)
			{


				//cout << " discover: " << u.value;                  /////////////////////////////////////////
				u.color = Gray;
				set_color(g, u, 1);
				stack.push_back(u);
				g.NumDiscover++;
				break; // if continue search all neighbors, became BFS 
			}


		}
		//cout << "\n";                                  /////////////////////////////////////////
	}

}

void DFS_test(Graph& g, int start) {
	for (int i = start; i < g.nVertexs; i++)
	{
		if (g.NumDiscover>=g.nVertexs)
		{
			cout << "g.NumDiscover: " << g.NumDiscover << endl;
			return;
		}
		Vertex curr = g.vertexs[i];
		//cout << curr.color;
		if (curr.color==White)
		{
			DFS_visit_test(g, curr);
		}
	}
	cout << "g.NumDiscover: " << g.NumDiscover<<endl;
}
void DFS(Graph& g, int start) {
	for (int i = start; i < g.nVertexs; i++)
	{

		Vertex curr = g.vertexs[i];
		//cout << curr.color;
		if (curr.color == White)
		{
			DFS_visit(g, curr);
		}

		if (g.NumDiscover >= g.nVertexs)
		{
			//cout << "g.NumDiscover: " << g.NumDiscover  << endl;
			return;
		}
	}
	//cout << "g.NumDiscover: " << g.NumDiscover << endl;
}

void DFS_optimize1(Graph& g, int threads) {
	int threadId = omp_get_thread_num();
	auto start = high_resolution_clock::now();
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);

	start = high_resolution_clock::now();

	for (int i = threadId; i < g.nVertexs; i=i+threads)
	{
		

		Vertex curr = g.vertexs[i];
		//cout << "\nat: " << curr.value<<endl;           //////////////////////////////////
		//printf("threadId: %d  ：  i is %d, at %d, color %d\n", omp_get_thread_num(), i, curr.value, curr.color);//////////////////////////////////
		if (curr.color == White)
		{
			//printf("threadId: %d  ：  i is %d, at %d, color %d\n", threadId, i, curr.value, curr.color);//////////////////////////////////
			//DFS_visit_test(g, curr);
			//DFS_visit(g, curr);
			DFS_visit_unOrder(g, curr, threadId, threads);
		}

		if (g.NumDiscover >= g.nVertexs)
		{
			//cout << "g.NumDiscover: " << g.NumDiscover  << endl;

			// stop = high_resolution_clock::now();
			// duration = duration_cast<milliseconds>(stop - start);
			// std::cout << threadId << " threadId spends (ms): " << duration.count() << endl;

			// printf("%d threadId spends (ms): %ld\n", threadId, duration.count());
			return;
		}
	}
	// stop = high_resolution_clock::now();
	// duration = duration_cast<milliseconds>(stop - start);
	// std::cout << threadId << " threadId spends (ms): " << duration.count() << endl;


	// printf("%d threadId spends (ms): %ld\n", threadId, duration.count());
	//cout << "g.NumDiscover: " << g.NumDiscover << endl;
}

//void DFS_optimize2_dynamic(Graph& g) {
//	int start = 0, end = g.nVertexs;
//#pragma omp parallel for schedule(dynamic)
//	for (int i = 0; i <end; i++)
//	{
//
//		Vertex curr = g.vertexs[i];
//		//cout << curr.color;
//		if (curr.color == White)
//		{
//			DFS_visit(g, curr);
//		}
//
//		if (g.NumDiscover >= g.nVertexs)
//		{
//			//cout << "g.NumDiscover: " << g.NumDiscover  << endl;
//			return;
//		}
//	}
//	//cout << "g.NumDiscover: " << g.NumDiscover << endl;
//}
//Experiment，run dfs 3 times pre_heat, and test  PDFS 10 time's averenge time
void Experiment(const char* filename, int threads) {
	auto start = high_resolution_clock::now();
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	int totalTime = 0;
	for (int i = 0; i < 13; i++)
	{
		Graph g = load_graph(filename);
		initialize_vertex(g);

		omp_set_num_threads(threads);
		start = high_resolution_clock::now();
		//////////////////////////////////////'
#pragma omp parallel 
		{
			
			//for (int i = 0; i < 2; i++)
			{
				//cout << "thread: " << omp_get_thread_num() << "---------------\n";
				DFS(g, omp_get_thread_num() * g.nVertexs / threads);

			}
		}

		//////////////////////////////////
		stop = high_resolution_clock::now();
		duration = duration_cast<milliseconds>(stop - start);
		

		if (i>2) {
			totalTime += duration.count();
			//printf("i %d   duration.count(): %d ms\n",i , duration.count());
			//cout<< duration.count()<<endl;
		}
	}

	cout << threads << " threads parallelized  time spends (ms): " << totalTime / 10 << endl;


//	totalTime = 0;
//	for (int i = 0; i < 13; i++)
//	{
//		Graph g = load_graph(filename);
//		initialize_vertex(g);
//
//		omp_set_num_threads(threads);
//		start = high_resolution_clock::now();
//		//////////////////////////////////////'
//#pragma omp parallel 
//		{
//
//			//for (int i = 0; i < 2; i++)
//			{
//				//cout << "thread: " << omp_get_thread_num() << "---------------\n";
//
//				DFS_optimize1(g, threads);
//
//			}
//		}
//
//		//////////////////////////////////
//		stop = high_resolution_clock::now();
//		duration = duration_cast<milliseconds>(stop - start);
//
//
//		if (i > 2) {
//			totalTime += duration.count();
//			//printf("i %d   duration.count(): %d ms\n",i , duration.count());
//			//cout<< duration.count()<<endl;
//		}
//	}
//
//	cout << threads << " DFS_optimize1(g, threads);  time spends (ms): " << totalTime / 10 << endl;
//


}







// int main() {
int main(int argc, char *argv[])  
{
    const char* filename = (argv[1]);
    int threads = atoi(argv[2]);
	cout<<"\n\n-----------------------------------------------------------------------------------------------\n";
	std::cout << filename<<" threads:"<<threads<<endl;
	//Graph g(10, true);
	//g.addEdge(0, 1);
	//g.addEdge(0, 4);
	//g.addEdge(1, 2);
	//g.addEdge(4, 2);
	//g.addEdge(4, 7);
	//g.addEdge(2, 7);
	//g.addEdge(2, 5);
	//g.addEdge(2, 3);
	//g.addEdge(3, 5);
	//g.addEdge(3, 6);
	//g.addEdge(5, 8);
	//g.addEdge(7, 8);
	//g.addEdge(8, 9);
	//cout << "graph g: \n";
	//initialize_vertex(g);
	//print_graph(g);


	// int threads = 56;
	// const char* filename= "com-lj.ungraph.txt";
	cout<<filename<<endl;
	Graph g= load_graph_cout(filename); //
	// email-Eu-core.txt  toy.dat test-1024.txt Wiki-Vote.txt test-5.5kEdges.txt   p2p-Gnutella31.txt 
	//soc-pokec-relationships.txt  com-orkut.ungraph.txt com-lj.ungraph.txt  Email-EuAll.txt

	initialize_vertex(g); // set all color as 0, sources as -1

#if 01           // for debug
	auto start = high_resolution_clock::now();
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);

	//print_graph(g);
	start = high_resolution_clock::now();
	//////////////////////////////////////'
	//DFS_test(g, 0);
	DFS(g,  0);
	//////////////////////////////////
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	std::cout << "\n baseline time spends (ms): " << duration.count() << endl;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	g = load_graph(filename); //
	// email-Eu-core.txt  toy.dat test-1024.txt Wiki-Vote.txt test-5.5kEdges.txt
	initialize_vertex(g); // set all color as 0, sources as -1
	omp_set_num_threads(threads);
	start = high_resolution_clock::now();
	//////////////////////////////////////'
#pragma omp parallel 
	{
			DFS(g, omp_get_thread_num()*g.nVertexs/ threads);
	}

	//////////////////////////////////
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	cout << "\n parallelized time spends (ms): " << duration.count() << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		 g = load_graph(filename);
		initialize_vertex(g);

		omp_set_num_threads(threads);
		start = high_resolution_clock::now();
		//////////////////////////////////////'
#pragma omp parallel 
		{

			//for (int i = 0; i < 2; i++)
			{
				//cout << "thread: " << omp_get_thread_num() << "---------------\n";

				DFS_optimize1(g, threads);

			}
		}

		//////////////////////////////////
		stop = high_resolution_clock::now();
		duration = duration_cast<milliseconds>(stop - start);
		cout << "\n parallelized DFS_optimize0 time spends (ms): " << duration.count() << endl;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	g = load_graph(filename); //
	// email-Eu-core.txt  toy.dat test-1024.txt Wiki-Vote.txt test-5.5kEdges.txt
	initialize_vertex(g); // set all color as 0, sources as -1
	omp_set_num_threads(threads);
	start = high_resolution_clock::now();
	//////////////////////////////////////'
#pragma omp parallel 
	{
			DFS(g, omp_get_thread_num()*g.nVertexs/ threads);
	}

	//////////////////////////////////
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	cout << "\n parallelized time spends (ms): " << duration.count() << endl;




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	g = load_graph_cout(filename);
	initialize_vertex(g);
	//print_graph(g);
	omp_set_num_threads(threads);
	start = high_resolution_clock::now();
	//////////////////////////////////////'
	int  end = g.nVertexs;
#pragma omp parallel for schedule(dynamic,50)   //////////////////////////////////////////////////////////////schedule guided  dynamic
	for (int i = 0; i < end; i++)
	{

		Vertex curr = g.vertexs[i];
		//cout << curr.color;
		if (curr.color == White)  // (curr.color == White)
		{
			DFS_visit(g, curr);
		}

		if (g.NumDiscover >= g.nVertexs)
		{
			//cout << "g.NumDiscover: " << g.NumDiscover  << endl;
			continue;
		}
	}
	//////////////////////////////////
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	std::cout << threads << " threads DFS_optimize1_dynamic(g); time spends (ms): " << duration.count() << endl;

	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	g = load_graph_cout(filename);
	initialize_vertex(g);
	//print_graph(g);
	omp_set_num_threads(threads);
	start = high_resolution_clock::now();
	//////////////////////////////////////'
	 end = g.nVertexs;
#pragma omp parallel for schedule(guided)   //////////////////////////////////////////////////////////////schedule guided  dynamic
	for (int i = 0; i < end; i++)
	{

		Vertex curr = g.vertexs[i];
		//cout << curr.color;
		if (curr.color == White)  // (curr.color == White)
		{
			DFS_visit(g, curr);
		}

		if (g.NumDiscover >= g.nVertexs)
		{
			//cout << "g.NumDiscover: " << g.NumDiscover  << endl;
			continue;
		}
	}
	//////////////////////////////////
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	std::cout << threads << " threads DFS_optimize2_guided(g); time spends (ms): " << duration.count() << endl;

#pragma omp barrier

	g = load_graph_cout(filename);
	initialize_vertex(g);
	//print_graph(g);
	omp_set_num_threads(threads);
	start = high_resolution_clock::now();

	int numOfthread = omp_get_num_threads();
	int* workload = new int[threads]();


#pragma omp parallel for schedule(guided)
	for (int i = 0; i < end; i++)
	{

		Vertex curr = g.vertexs[i];
		//cout << curr.color;
		if (curr.color == White)  // (curr.color == White)
		{
			//DFS_visit(g, curr);
			DFS_visit_workload(g, curr, workload);
		}

		if (g.NumDiscover >= g.nVertexs)
		{
			//cout << "g.NumDiscover: " << g.NumDiscover  << endl;
			//for (int i = 0; i < sizeof(workload)/sizeof(workload[0]); i++)
			//{
			//	cout << i << " thread's workload: " << workload[i] << endl;
			//}
			continue;
		}
	}
#pragma omp barrier
	//////////////////////////////////
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	std::cout << threads << " threads DFS_optimize2_dynamic(g); time spends (ms): " << duration.count() << endl;
	for (int i = 0; i < threads;  i++)
	{
		cout << i << " thread's workload: " << workload[i] << endl;
	}




#else  // use openmp
	auto start = high_resolution_clock::now();
	auto stop = high_resolution_clock::now();


	auto duration = duration_cast<milliseconds>(stop - start);
	start = high_resolution_clock::now();
	//////////////////////////////////////'
	DFS(g, 0);
	//////////////////////////////////
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	std::cout << "\n serial time spends (ms): " << duration.count() << endl;

	g = load_graph(filename); //
	// email-Eu-core.txt  toy.dat test-1024.txt Wiki-Vote.txt test-5.5kEdges.txt
	initialize_vertex(g); // set all color as 0, sources as -1

	//MyVisitor vis;
	
	omp_set_num_threads(threads);
	cout << "\nomp_get_thread_num: " << omp_get_thread_num() << "--------------------------------\n" << endl;
	start = high_resolution_clock::now();
	//////////////////////////////////////'
#pragma omp parallel 
	{
		//for (int i = 0; i < 2; i++)
		{

			//cout << "thread: " << omp_get_thread_num() << "---------------\n";
			//printf("thread: %d start at %d\n", omp_get_thread_num(), omp_get_thread_num() * g.nVertexs / 2);
			DFS(g, omp_get_thread_num()*g.nVertexs/ threads);
			
		}
	}

	//////////////////////////////////
	stop = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(stop - start);
	cout << "\n parallelized  time spends (ms): " << duration.count() << endl;


	cout << "Experiment "  << "----------------------------------------------------------------------------------------\n" << endl;
	for (int i = 1; i < 7; i++)
	{

		Experiment(filename, i);
	}

#endif // 0
	//g.addEdge(6, 7);
	////vector<int> bv = g.BFS(0);
	//cout << "BFS: ";
	//for (int j = 0; j < bv.size(); j++)
	//	cout << bv[j] << " ";

	//cout << endl;
	//cout << "DFS: ";
	//Graph g1(6, false);
	//g1.addEdge(0, 1);
	//g1.addEdge(0, 4);
	//g1.addEdge(0, 5);
	//g1.addEdge(1, 5);
	//g1.addEdge(4, 5);
	//g1.addEdge(5, 2);
	//g1.addEdge(5, 3);
	//g1.addEdge(2, 3);
	//vector<int> route = g1.DFS(0);
	//for (int i = 0; i < route.size(); i++)
	//	cout << route[i] << " ";
	//cout << endl;
	//char ch;
	//cin >> ch;
	return 0;
}


