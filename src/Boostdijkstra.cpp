#include "Boostdijkstra.hpp"

#include <stdexcept>
#include <iostream> 
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each

#include <boost/graph/adjacency_list.hpp> 
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>     
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _MPI
#include <mpi.h>
#endif

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

using namespace std; 
using namespace boost;



//Foo::Foo(int _a, int _b): a(_a), b(_b){
//    cout << "C++ side, constructor" << endl;
//}
//
//Foo::~Foo(){
//    cout << "C++ side, destructor" << endl;
//}
//
//int Foo::bar(int c) const{
//    return a + c;
//}
//
//double Foo::baz(double d) const{
//    return d + b;
//}

//void foo_speaker(string s){
////    Foo f(4, 2);
////    cout << s << " Foo(4, 2).bar(3) is: " <<  f.bar(3) << endl;
//     cout << s << endl;
//}

//
////void dijkstra(string s) 
//void dijkstra_cpp_mpi(int* NODES_LENGTH, int* GRID_OUT_LENGTH, int* cur_edge_node_array1, int* cur_edge_node_array2, double* edge_dist_array, int* NODES2CAL_LENGTH, int* nodes2cal_array, double* coord_ISOMAP) 
//{
//    cout << "Begin distributed dijkstra" << endl;
//    ierr = MPI_Init ( &argc, &argv );
//
//
//}
//
void dijkstra_cpp(int* NODES_LENGTH, int* GRID_OUT_LENGTH, int* cur_edge_node_array1, int* cur_edge_node_array2, double* edge_dist_array, int* NODES2CAL_LENGTH, int* nodes2cal_array, double* coord_ISOMAP) 
{ 
    cout << " [CPP] Begin dijkstra" << endl;

//    cout << *NODES_LENGTH <<" " << *GRID_OUT_LENGTH << " " << *(cur_edge_node_array1+0) << " " <<  *(cur_edge_node_array2+0) << " " << *(edge_dist_array) << endl;
    //cout << *(coord_ISOMAP+0) << " "<< *(coord_ISOMAP+1) << " " << *(coord_ISOMAP+*NODES_LENGTH) << endl;

//     start the timer
    time_t start_time, end_time;
    time_t diff_time;


    //typedef adjacency_list < listS, vecS, directedS, no_property, property < edge_weight_t, float > > graph_t;
    typedef adjacency_list < listS, vecS, undirectedS, no_property, property < edge_weight_t, float > > graph_t;
    typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
    typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
    typedef std::pair<int, int> Edge;

    long int num_nodes;
    long int num_edges;
    int nodes2cal;
    char gfile[20] ; //graph data file
    char nfile[20] ; //nodes to cal dist for data file
    char  buffer[128] ;
    istringstream instream ;

   


    //string gfile_str = "gird 1.6m nodes.out";
    //string gfile_str = "girdsmall.out";
    //string gfile_str = "grid.out";
    string gfile_str = "grid.out";
    strcpy(gfile,gfile_str.c_str());

//    ifstream infile( gfile, ios::in );
//    infile.getline(buffer, 128);
//    instream.clear() ;
//    instream.str(buffer) ;
//    instream >> num_nodes >> num_edges;
    num_nodes = *NODES_LENGTH;
    num_edges = *GRID_OUT_LENGTH;
    nodes2cal = *NODES2CAL_LENGTH;


    Edge* edge_array;   // Pointer to edge, initialize to nothing.
    edge_array = new Edge[num_edges];  // Allocate n ints and save ptr in a.
    float* weights;   // Pointer to edge, initialize to nothing.
    weights = new float[num_edges];  // Allocate n ints and save ptr in a.
    int num_arcs = num_edges;
    int n1,n2;

    //for(int k=0;k<num_nodes*nodes2cal;k++)
    //	coord_ISOMAP[k]=0.0;


    cout << " [CPP] Loading graph";  cout << endl;cout << endl;  cout << endl;
    cout << " [CPP] Found ";cout << num_nodes; cout << " nodes, and "; cout << num_edges; cout << " edges";cout  << endl; 

    for (int i=0; i<num_edges; ++i) 
    { 
        //get a line from the file
        //if (num_edges/2 == i || (num_edges+1)/2 == i )
        //{
        //    cout << "half way done reading the graph ... please wait a bit longer ...";  cout << endl; 
        //}

//        infile.getline(buffer, 128);
//        instream.clear() ;
//        instream.str(buffer) ;
//        instream >> n1 >> n2 >> weights[i];
        n1 = *(cur_edge_node_array1 + i);
        n2 = *(cur_edge_node_array2 + i);
        weights[i] = (float)(*(edge_dist_array + i));
        edge_array[i]=Edge(n1-1,n2-1);
    }

//     _strdate( dateStr1);
//     _strtime( timeStr1 );

    cout << " [CPP] Graph read in"; cout << endl;  cout << endl; 
    //cout << " [CPP] Setting up Dijkstra call for the graph"; cout << endl; 

    graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);

    //cout << "graph_t g created"; cout << endl; 

    property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);

    //cout << "property_map weightmap created"; cout << endl; 

    std::vector<vertex_descriptor> p(num_vertices(g));

    //cout << "vector<vertex_descriptor> p created"; cout << endl; 

    std::vector<float> d(num_vertices(g));

    //cout << "vector<float> d created"; cout << endl; 
    
    property_map<graph_t, vertex_index_t>::type indexmap = get(vertex_index, g);

    //cout << "property_map indexmap created"; cout << endl; 



//need to find out what distances are required?
//in file 'node2cal.out'
//first line is number of nodes to cal
//all other lines contain the index of the node to calculate
//remember in c++ the ind = indFORTRAN-1

    int cur_node;
    //string nfile_str = "nodes2cal.out";
    //strcpy(nfile, nfile_str.c_str());
    //ifstream infile2( nfile, ios::in ); 
    //infile2.getline(buffer, 128);
    //instream.clear() ;
    //instream.str(buffer) ;
    //instream >> nodes2cal;

    //nodes2cal = *NODES2CAL_LENGTH;
    cout << " [CPP] Working out distances to ";cout << nodes2cal; cout << " nodes"; cout << endl; 
    //cout << "NODE:"; cout << endl; 


#ifndef _OPENMP
    //now need to be able to write out the distances
//    ofstream outfile;
//    outfile.open ("dist_cpp.out");
    diff_time=0;
    for (int i=0; i<nodes2cal; ++i) 
    //for (int i=0; i<4; ++i) 
    {
//        cout << i + 1 << flush; 
//        cout << endl; 

        //infile2.getline(buffer, 128);
        //instream.clear() ;
        //instream.str(buffer) ;
        //instream >> cur_node;
        cur_node=nodes2cal_array[i];

        vertex_descriptor s = vertex(cur_node-1, g);  //node number output from fortran program will always be +1 (C++ starts at 0 not 1)
         start_time = time(NULL);

//        dijkstra_shortest_paths(g, s, &p[0], &d[0], weightmap, indexmap, 
//                          std::less<float>(), closed_plus<float>(), 
//                          (std::numeric_limits<float>::max)(), 0,
//                          default_dijkstra_visitor());


	  	dijkstra_shortest_paths(g, s,  
	       		predecessor_map(
	        	 	make_iterator_property_map(p.begin(), get(vertex_index, g))).
	       		distance_map(
	        		make_iterator_property_map(d.begin(), get(vertex_index, g))).
	  		weight_map(get(edge_weight, g)).
	  		vertex_index_map(get(vertex_index, g)).
			distance_compare(std::less<float>()).
			distance_combine(closed_plus<float>()).
			distance_inf((std::numeric_limits<float>::max)()).
			distance_zero(0).
			visitor(default_dijkstra_visitor())
		);



        //write these distances out
	// Descomentar para escribir resultado
	// Comentar para medir speedup y escalabilidad
	//int inodes_length=i*(*NODES_LENGTH);
	long int inodes_length=i * num_nodes;
        for (int j=0; j<num_nodes; ++j) 
        {
            //outfile << d[j] <<  endl ; 
            //*(coord_ISOMAP+inodes_length+j)=d[j];
            coord_ISOMAP[inodes_length+j]=d[j];
        }
	

         end_time=time(NULL);
         diff_time=diff_time + end_time-start_time;
    }

    cout << " [CPP] Time for one path:  " << float(diff_time)/float(nodes2cal) <<  " s" << endl ;
    cout << " [CPP] Done dijkstra"; 
    //outfile.close();

#else
    //now need to be able to write out the distances

    int num_threads_tmp=0;
    int num_threads=0;
    int this_thread=0;
    //cout << omp_get_max_threads() << endl;;
    //omp_set_num_threads(16);
#pragma omp parallel default(none) firstprivate(this_thread) shared(num_threads,cout)
{
    //this_thread = omp_get_thread_num();
    //num_threads = omp_get_num_threads();
    #pragma omp critical
    {
    	this_thread=num_threads;
        //cout << " [CPP][OMP] this_thread=" << this_thread << " num_threads=" << num_threads << endl;
    	num_threads++;
    }
    #pragma omp barrier
}

    num_threads_tmp=num_threads;
    cout << " [CPP][OMP] num_threads=" << num_threads_tmp << endl; 

    int cur_node_array[nodes2cal];
    //string stringout = "dist_cpp.out"; 

    for (int i=0; i<nodes2cal; ++i) 
    {
        //infile2.getline(buffer, 128);
        //instream.clear() ;
        //instream.str(buffer) ;
        //instream >> cur_node_array[i];

        cur_node_array[i] = nodes2cal_array[i];

    }

    diff_time=0;
    start_time = time(NULL);

    this_thread=0;
    num_threads=0;
    int ichunk;
    int imin;
    int imax;
    int i,j;
    long int inodes_length;

#pragma omp parallel default(none) firstprivate(this_thread,cur_node,p,d,ichunk,imin,imax,i,j,inodes_length) \
                                   shared(num_nodes,num_threads,num_threads_tmp,nodes2cal,cout,cur_node_array,g,weightmap,indexmap,coord_ISOMAP,NODES_LENGTH)
{

    
    //this_thread = omp_get_thread_num();
    //num_threads = omp_get_num_threads();
    #pragma omp critical
    {
    	this_thread=num_threads;
        cout << " [CPP][OMP] this_thread=" << this_thread << " num_threads=" << num_threads << endl;
        //cout << this_thread << " " << num_threads << endl;
    	num_threads++;
    }
    
    #pragma omp barrier
    //while(num_threads!=num_threads_tmp);

    char num2str[21];

//    if(num_threads>1000){
//        cout << "Too many threads (max=1000). Exit" << endl;
//        return 1;
//    }

    //if(this_thread<10)
    //	sprintf(num2str,"_00%d",this_thread);
    //else if(this_thread<100)
    //	sprintf(num2str,"_0%d",this_thread);
    //else if(this_thread<1000)
    //	sprintf(num2str,"_%d",this_thread);

    //string stringoutLocal = stringout + string(num2str);

    //cout << this_thread << ":" << stringoutLocal << endl;

//    ofstream outfile;
//    outfile.open(stringoutLocal.c_str());

//#pragma omp for nowait
    ichunk=int((nodes2cal-num_threads_tmp+1)/num_threads_tmp)+1; 
    imin=this_thread*ichunk;
    imax=MIN((this_thread+1)*ichunk,nodes2cal);
    if(this_thread==num_threads_tmp-1){
    	imax=nodes2cal;
    }

//    #pragma omp critical
//    {
//    cout << " [CPP][OMP] " << this_thread << ": ichunk=" << ichunk << " imin=" << imin << " imax=" << imax << endl;
//    }
    #pragma omp barrier

    //for (int i=0; i<nodes2cal; ++i) 
    for (i=imin; i<imax; ++i) 
    {
        if(this_thread==0){cout << " [CPP][OMP] Estimated progress: " << i << "/" << ichunk << endl;}

        cur_node=cur_node_array[i];

        vertex_descriptor s = vertex(cur_node-1, g);  //node number output from fortran program will always be +1 (C++ starts at 0 not 1)

//      	dijkstra_shortest_paths(g, s,  
//           		predecessor_map(
//            	 	make_iterator_property_map(p.begin(), get(vertex_index, g))).
//           		distance_map(
//            		make_iterator_property_map(d.begin(), get(vertex_index, g))).
//      			weight_map(get(edge_weight, g)).
//      			vertex_index_map(get(vertex_index, g)).
//    			distance_compare(std::less<float>()).
//    			distance_combine(closed_plus<float>()).
//    			distance_inf((std::numeric_limits<float>::max)()).
//    			distance_zero(0).
//    			visitor(default_dijkstra_visitor())
//    	);
      	dijkstra_shortest_paths(g, s,  
           		predecessor_map(
            	 	make_iterator_property_map(p.begin(), indexmap)).
           		distance_map(
            		make_iterator_property_map(d.begin(), indexmap)).
      			weight_map(weightmap).
      			vertex_index_map(indexmap).
    			distance_compare(std::less<float>()).
    			distance_combine(closed_plus<float>()).
    			distance_inf((std::numeric_limits<float>::max)()).
    			distance_zero(0).
    			visitor(default_dijkstra_visitor())
    	);




        //write these distances out
        // Descomentar para escribir resultado
        // Comentar para medir speedup y escalabilidad
        
        //for (int j=0; j<num_nodes; ++j) 
        //{
        //    outfile << d[j] <<  endl ; 
        //}  
    	//inodes_length=i*(*NODES_LENGTH);
	inodes_length=i * num_nodes;
        //cout << " [CPP][OMP] this_thread="<< this_thread <<" inodes_length=" << inodes_length << " num_nodes="<< num_nodes << endl;
	for (j=0; j<num_nodes; j++) 
        {
            //outfile << d[j] <<  endl ; 
            //*(coord_ISOMAP+inodes_length+j)=d[j];
            coord_ISOMAP[inodes_length+j]=d[j];
        }    

    }


//    outfile.close();

}

    end_time=time(NULL);
    diff_time=diff_time + end_time-start_time;
    cout << " [CPP] Time for one path:  " << float(diff_time)/float(nodes2cal) <<  " s" << endl ;
    cout << " [CPP] Done dijkstra" << endl; 


#endif

    delete(edge_array) ;
    delete(weights) ;

//    return 0; 
} 

