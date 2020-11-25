
#include "trail.h"
#include <limits>
using namespace std;

void Trail::Run()
{
	printf("begin trail, threshold %lf, number of file %d\n",e,n);
	int *trail_set =  (int *)malloc(sizeof(int) * n);
	memset(trail_set,0,n);
	for(int i = 0; i < n; i++)
		trail_set[i] = 0;

	float prev_min = dist_matrix[0][1];
	int start = 0, current = 1;
	for(int i = 0; i < n-1; i++)
	{
		for(int j = i+1; j < n; j++)
		{
			if( dist_matrix[i][j] < prev_min )
			{
				start = i;
				current = j;
				prev_min = dist_matrix[i][j];  
			}
		}
	}

	trail_set[start] = 1; 
	trail_set[current] = 1;
	// int start_max = FindMin(start,trail_set);
	// int cur_max = FindMin(current,trail_set);
	// //printf("test: %.4f %.4f\n",w_matrix[start][start_max],w_matrix[current][cur_max]);
	// if ( dist_matrix[start][start_max] < dist_matrix[current][cur_max] )
	// {
	// 	int buffer = start;
	// 	start = current;
	// 	current = buffer;
	// }

	group *cur_group = (group *)malloc(sizeof(group));
	kv_init(*cur_group);
	kv_push(int, *cur_group, start);
	kv_push(int, *cur_group, current);

	int count = 2;
	vector <int> test_result;
	test_result.push_back(start);
	test_result.push_back(current);
	float average_dist = 0.0;
	for (int i = 0; i < n-1; i++)
	{
		for(int j = i+1; j < n; j++)
			average_dist += dist_matrix[i][j];
	}

	//printf("test dist %.4f, test trail %d\n",dist_matrix[60][82],trail_set[82]);

	dist_threshold = 0.5*average_dist/(n*(n-1));
	
	printf("distance: %d,%d,%lf\n",start,current,dist_matrix[start][current]);
	while (count < n)
	{
		int next = FindMin(current,trail_set);
		printf("distance: %d,%d,%lf\n",current,next,dist_matrix[current][next]);
		trail_set[next] = 1;
		if( dist_matrix[current][next] <= e * (0.9+2/(cur_group->n+9))* prev_min || 
			(dist_matrix[current][next] < 1.2*dist_threshold && cur_group->n <=4) )
		{
			kv_push(int, *cur_group, next);
		}
		else
		{
			groups.push_back(cur_group);
			cur_group = (group *)malloc(sizeof(group));
			kv_init(*cur_group);
			kv_push(int, *cur_group, next);
		}

		test_result.push_back(next);
		prev_min = dist_matrix[current][next];
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

int Trail::FindMin(int id, int *trail_set)
{
	float min = numeric_limits<float>::max();
	int result = -1;
	if(id == 60)
	{

		for(int i = 0 ; i < n; i++)
			printf("%d ",trail_set[i]);
		printf("\n");
		for(int i = 0 ; i < n; i++)
		{

			if( trail_set[i] != 1 )
				printf("remaining point:%d,%.4lf",i,dist_matrix[id][i]);
		}
	}
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
		float prev_dist = numeric_limits<float>::max();
		int count = 0, cur_num = 0;
		while( count < (groups[i]->n)-2 && groups[i]->n > 3)
		{
			//printf("step %d ",count);
			int cur = groups[i]->a[count], next = groups[i]->a[count+1];
			float cur_dist = dist_matrix[cur][next];
			double threshold = 1 + 0.25*(4-cur_num); 
			if(cur_num <= 2 && 	(cur_dist <= threshold * prev_dist || cur_dist <1.2*dist_threshold))
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
			prev_dist = dist_matrix[cur][next];
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