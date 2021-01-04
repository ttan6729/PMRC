#include <bits/stdc++.h> 

#include <bitset>  
#include <assert.h>
#include <stdint.h>
#include <unistd.h>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <cmath>  
#include <assert.h>
#include <omp.h>
#include "kvec.h"
//#include <vectors>
// A structure to represent a subset for union-find 
typedef struct
{
    int parent; //root of current vertic
    int rank;  //rank of current vertic
}   subset;

typedef struct   
{  
    int vertices[2]; //first is id of src, second is id of destinaiton
    double weight;  
} Edge;

typedef struct {	int n,m;	Edge *a;	} Edges; //vector of edge

int EdgeComp(const void* a, const void* b);
int check_subsets(subset *subsets, int i);

class Graph  
{  
public: 
    int V;   
    Edges *es;
    double threshold=1.5; 
    Edges *result;
    Graph(int _V, double _threshold)
    {
    	V = _V;
        threshold = _threshold;
    	es = (Edges *)malloc(sizeof(Edges));
    	kv_init(*es);
    }
    void AddEdge(Edge e)	{	kv_push(Edge,*es,e);	}
    void Union(subset *subsets, int x, int y);  
    void MST();
    void Split(char *output_fp);
    ~Graph()
    {
        kv_destroy(*es);

    }
};  


class PrimMST
{
public:
    float **dist_matrix; //
    int V; //number of vertices
    int *parent; //constructed MST
    PrimMST(float **_matrix, int _V): dist_matrix(_matrix), V(_V) 
    {
        parent = (int *)malloc(sizeof(int) * _V);
    }
    
    void RunMST();
    ~PrimMST()
    {

    }
};