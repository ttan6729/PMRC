#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>
#include "creads.h"
//#include "config.h"
#include "cseq.h"
#include "mst.h"
#include "trail.h"
using namespace std;

void RunMST(float **dist_matrix);
float **process(char *dir, int k, char *matrix_fp);
float **ReadMatrix(char *file_name);

char *file_list_path;
int n_threads=12;
double threshold=1.2;
int selected_number=100;
char *output_fp;
vector<string> file_list;

int main (int argc, char **argv)
{	
	int _k = 7;
	int c;
	int if_read = 0;
	char *matrix_fp=NULL;
	while ((c = getopt(argc, argv, "r:o:t:k:e:s:m:z:")) != -1)	
    switch (c)
    {
      	case 'r':
        file_list_path = optarg;
        break;
     	case 't':
        n_threads = atoi(optarg);
        break;
        case 'o':	
     	output_fp = optarg;
     	break;
      	case 'k':
        _k = atoi(optarg);
        break;
        case 's':
        selected_number = atoi(optarg);
        break;
      	case 'e':
      	threshold = atof(optarg);
      	break;
      	case 'm':
      	file_list_path = optarg;
      	if_read = 1;
      	case 'z':
      	matrix_fp = optarg;
      	break;      	
    }
    char output_name[1024];
    sprintf(output_name,"%s/cluster",output_fp);
    if( if_read )
    {
    	float **dist_matrix = ReadMatrix(file_list_path);
    	printf("file number:%d\n",file_list.size());
    	Trail *trail = new Trail(file_list.size(),dist_matrix);
    	trail->Run();
    	trail->WriteResult(output_fp);
    	//RunMST(dist_matrix);
    	return 0;
    }
    else 
    {	
    	float **dist_matrix = process(output_fp,_k,matrix_fp); 
    	Trail *trail = new Trail(file_list.size(),dist_matrix);
    	trail->Run();
    	trail->WriteResult(output_name);
    }
	return 0;
}

float **process(char *dir,int k, char *matrix_fp)
{
	printf("input file list path: %s\n",file_list_path);
	fstream file;
	file.open(file_list_path,ios::in); 
   	int n = 0;
   	if (file.is_open())
   	{   //checking whether the file is open
      string line;
      while(getline(file, line))
      { 
        file_list.push_back(line);
        n++;
      }
      file.close(); //close the file object.
    }
    else
    {
    	printf("failed to open %s\n",file_list_path);
    	exit(1);
    }
	float **vectors = pre_process(file_list,output_fp,k,selected_number,matrix_fp);
	return vectors;
}

void RunMST(float **dist_matrix)
{
	printf("bgein minimum tree\n");
	// Edge e;
	// Graph* graph = new Graph(file_list.size(),threshold);
	// for(int i = 0; i < file_list.size()-1; i++)
	// {
	// 	for(int j = i+1; j < file_list.size(); j++)
	// 	{
	// 		e.vertices[0] = i;
	// 		e.vertices[1] = j;
	// 		e.weight = dist_matrix[i][j];
	// 		graph->AddEdge(e);
	// 	}
	// }
	PrimMST *pMST = new PrimMST(dist_matrix,file_list.size());
	printf("begin MST\n");
	pMST->RunMST();
	//graph->Split(output_fp);
}

float **ReadMatrix( char *file_name)
{

	ifstream f(file_name);
	string line,value;
	int n =	count(istreambuf_iterator<char>(f), istreambuf_iterator<char>(),'\n')-1;
	f.close();
	printf("read file %s, line number:%d\n",file_name,n);
	// int n = 16, m = 4096;
	float **vectors = (float **)malloc(sizeof(float *)*n);
	for(int i = 0; i < n; i++)
	{
		file_list.push_back("a");
		vectors[i] =  (float *)malloc(sizeof(float)*selected_number);
	}

 	ifstream myfile;
    myfile.open(file_name);
	getline(myfile,line);
	for(int i = 0; i < n; i++)
	{
		getline(myfile,line);
		stringstream linestream(line);
		getline(linestream,value,',');
		for(int j = 0; j < selected_number; j++)
		{
			getline(linestream,value,',');
			vectors[i][j] = stof(value);
		}
	}

	// for(int i = 0; i < n; i++)
	// {
	// 	float sum = 0.0;
	// 	for(int j = 0; j < selected_number; j++)
	// 		sum += vectors[i][j];
	// 	for(int j = 0; j < selected_number; j++)
	// 		vectors[i][j] = 100*vectors[i][j]/sum;
	// }
	myfile.close();	
	return vectors;
}
// void RunMST()
// {
// 	//test case
//     Graph* graph = new Graph(4,1.2);
//         // add edge 0-1  
//    	Edge e;
// 	e.vertices[0] = 0;  
//     e.vertices[1] = 1;  
//     e.weight = 10.0;  
//     graph->AddEdge(e);
  
// 	e.vertices[0] = 0;  
//     e.vertices[1] = 2;  
//     e.weight = 6.0;  
//     graph->AddEdge(e);

// 	e.vertices[0] = 0;  
//     e.vertices[1] = 3;  
//     e.weight = 5.0;  
//     graph->AddEdge(e);

// 	e.vertices[0] = 1;  
//     e.vertices[1] = 3;  
//     e.weight = 15.0;  
//     graph->AddEdge(e);  

// 	e.vertices[0] = 2;  
//     e.vertices[1] = 3;  
//     e.weight = 4.0; 
//     graph->AddEdge(e);

//     graph->MST();
//     graph->Split(output_fp);

//     delete graph;
//     return;
// }

