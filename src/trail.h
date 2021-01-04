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

typedef struct 
{
	int n;
	int m;
	int *a; //id of points
} group;

class Trail
{
public:
	char *file_name;
	std::vector<group *> groups;
 	int n; //number of points
 	float **dist_matrix;
 	float **w_matrix;
 	float e = 1.3;
 	float w_threshold, dist_threshold;
 	void Run();
 	int FindMax(int id, int *trail_set);
 	int FindMin(int id, int *trail_set);
 	void WriteResult(char *file_name);
 	void Split(group *g);
 	Trail(int _n, float **_dist_matrix, float _e = 1.3) 
 	{
 		e = _e;
 		n = _n;
 		dist_matrix = _dist_matrix;
 		// w_matrix = (float **)malloc(sizeof(float *) *n);
 		// w_matrix[0] = (float *)malloc(sizeof(float) *n*n);
 		// for(int i = 1; i < n; i++)
 		// 	w_matrix[i] = w_matrix[0] + n * i;
 		// for(int i = 0; i < n; i++)
 		// {
 		// 	w_matrix[i][i] = 1.0;
 		// 	for(int j = i+1; j < n; j++)
 		// 	{
 		// 		//printf("%.4f %.4f\n",dist_matrix[i][j],exp(-0.1*dist_matrix[i][j]));
 		// 		w_matrix[i][j] = 0.01+exp(-0.01*dist_matrix[i][j]);
 		// 		w_matrix[j][i] = w_matrix[i][j];
 		// 	}
 		// }
 		printf("dimension %d\n",n);

 	} 	

 	~Trail()
 	{
 		//for(int i = 0; i < n; i++)
 		//	free(w_matrix[i]);
 		//free(w_matrix);

 		for(int i = 0; i < n; i++)
 			free(dist_matrix[i]);
 		free(dist_matrix);
 		for(int i = 0; i < groups.size(); i++)
 			free(groups[i]);
 	}

};

