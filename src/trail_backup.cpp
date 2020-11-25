#include "trail.h"
#include <limits>
using namespace std;

void Trail::Run()
{
	// printf("test weight matrix\n");
	// for(int i = 0; i < n; i++)
	// {
	// 	for(int j = 0; j < n; j++)
	// 		printf("%.3f ",w_matrix[i][j]);
	// 	printf("\n");
	// }
	printf("begin trail, threshold %lf, number of files %d\n",e,n);
	int *trail_set =  (int *)malloc(sizeof(int) * n);
	memset(trail_set,0,n);

	float prev_max = 0.0;
	int start = 0, current = 1;
	for(int i = 0; i < n-1; i++)
	{
		for(int j = i+1; j < n; j++)
		{
			if(  w_matrix[i][j] > prev_max )
			{
				start = i;
				current = j;
				prev_max = w_matrix[i][j];  
			}
		}
	}

	trail_set[start] = 1; trail_set[current] = 1;
	int start_max = FindMax(start,trail_set);
	int cur_max = FindMax(current,trail_set);
	//printf("test: %.4f %.4f\n",w_matrix[start][start_max],w_matrix[current][cur_max]);
	if ( w_matrix[start][start_max] > w_matrix[current][cur_max] )
	{
		int buffer = start;
		start = current;
		current = buffer;
	}

	group *cur_group = (group *)malloc(sizeof(group));
	kv_init(*cur_group);
	kv_push(int, *cur_group, start);
	kv_push(int, *cur_group, current);

	int count = 2;
	vector <int> test_result;
	test_result.push_back(start);
	test_result.push_back(current);

	float average_weight = 0.0;
	for (int i = 0; i < n-1; i++)
	{
		for(int j = i+1; j < n; j++)
			average_weight += w_matrix[i][j];
	}

	w_threshold = 0.5*average_weight/(n*(n-1));
	printf("weight: %d,%d,%lf\n",start,current,w_matrix[start][current]);
	while (count < n)
	{
		int next = FindMax(current,trail_set);
		//printf("weight: %d,%d,%lf\n",current,next,w_matrix[current][next]);
		trail_set[next] = 1;
		if( w_matrix[current][next] >= e * (0.9-3/(cur_group->n+9))* prev_max || 
			(w_matrix[current][next] > w_threshold && cur_group->n <=4) )
		{
			//printf("a\n");
			kv_push(int, *cur_group, next);
		}
		else
		{
			//printf("b\n");
			groups.push_back(cur_group);
			cur_group = (group *)malloc(sizeof(group));
			kv_init(*cur_group);
			kv_push(int, *cur_group, next);
		}
		test_result.push_back(next);
		prev_max = w_matrix[current][next];
		current = next;
		count++;
	}
	groups.push_back(cur_group);
	printf("test result:\n");
	for(int i = 0; i < n; i++)
		printf("%d ", test_result[i]);
	printf("\n");
	free(trail_set);
	return;
}

int Trail::FindMax(int id, int *trail_set)
{
	float max = 0.0;
	int result = -1;
	for(int i = 0 ; i < n; i++)
	{
		if( w_matrix[id][i] > max)
		{
			if( i != id && trail_set[i] != 1 )
			{
				result = i;
				max = dist_matrix[id][i];
			}
		}
	}

	return result;
}


int Trail::FindMin(int id, int *trail_set)
{
	float min = numeric_limits<float>::max();
	int result = -1;
	for(int i = 0 ; i < n; i++)
	{
		if( dist_matrix[id][i] < min)
		{
			if( i != id && trail_set[i] != 1 )
			{
				result = i;
				min = dist_matrix[id][i];
			}
		}
	}

	return result;
}

void Trail::WriteResult(char *file_name)
{
	printf("result write to %s\n",file_name);
	ofstream fp;
	fp.open(file_name);
	printf("begin split, number of groups: %d\n", groups.size());
	int total = 0;
	for(int i = 0; i < groups.size(); i++)
	{
		total += groups[i]->n;
		printf("start %d, group size:%d\n",groups[i]->a[0],groups[i]->n);
		float prev_weight = 0.0;
		int count = 0, cur_num = 0;
		while( count < (groups[i]->n)-2 && groups[i]->n > 3)
		{
			//printf("step %d ",count);
			int cur = groups[i]->a[count], next = groups[i]->a[count+1];
			float cur_weight = w_matrix[cur][next];
			double threshold = 5/(10-cur_num); 
			if(cur_num <= 2 && 	(cur_weight >= threshold * prev_weight || cur_weight > w_threshold))
			{
				//printf("a, count %d\n",count);
			 	fp << cur << " ";
			 	cur_num++;
			}
		 	else
			{
				//printf("b, count %d\n",count);
			 	fp <<  "\n" << cur << " ";
			 	cur_num = 1;
			}
			count++;
			prev_weight = w_matrix[cur][next];
			if(cur_num >=2 && count >= ((groups[i]->n)-2) )
				fp << "\n";
		}
		for(int j = count; j < groups[i]->n; j++ )
			fp << groups[i]->a[j] << " ";
		fp << "\n"; 
	}
	printf("included file %d\n",total);
	fp.close();
	return;

}

void Trail::Split(group *g)
{

}
