#include "mst.h"
using namespace std;

int EdgeComp(const void* a, const void* b)
{  
    Edge* a1 = (Edge*)a;  
    Edge* b1 = (Edge*)b;  
    return a1->weight > b1->weight;  
}  

int CheckSubsets(subset *subsets, int i)
{  
    // find root and make root as parent of i   
    if (subsets[i].parent != i)  
        subsets[i].parent = CheckSubsets(subsets, subsets[i].parent);  
    return subsets[i].parent;  
}  

void Graph::Union(subset *subsets, int x, int y)  
{  
    int xroot = CheckSubsets(subsets, x);  
    int yroot = CheckSubsets(subsets, y);  
  
    // Attach smaller rank tree under root of high  
    // rank tree (Union by Rank)  
    if (subsets[xroot].rank < subsets[yroot].rank)  
        subsets[xroot].parent = yroot;  
    else if (subsets[xroot].rank > subsets[yroot].rank)  
        subsets[yroot].parent = xroot;  
  
    // If ranks are same, then make one as root and  
    // increment its rank by one  
    else
    {  
        subsets[yroot].parent = xroot;  
        subsets[xroot].rank++;  
    }  
}  

//construct MST using Kruskal's algorithm  
void Graph::MST()  
{ 
	if( V < 2)
	{
		printf("error not, enough vertices for minimum spanning tree\n");
		return;
	} 

    result =  (Edges *)malloc(sizeof(Edges)); // Tnis will store the resultant MST  
    kv_init(*result);
    int iter = 0; // An index variable, used for sorted edges  
  	
  	qsort(es->a, es->n, sizeof(Edge), EdgeComp);
  
    // Allocate memory for creating V ssubsets  
    subset *subsets = new subset[( V * sizeof(subset) )];  
  
    // Create V subsets with single elements  
    for (int v = 0; v < V; ++v)  
    {  
        subsets[v].parent = v;  
        subsets[v].rank = 0;  
    }  
  
    // Number of edges to be taken is equal to V-1  
    while ( result->n < V - 1 && iter < es->n)  
    {  
        Edge next_edge = es->a[iter++];  
  
        int x = CheckSubsets(subsets, next_edge.vertices[0]);  
        int y = CheckSubsets(subsets, next_edge.vertices[1]);  
  
        // If including this edge does't cause cycle,  
        // include it in result and increment the index  
        // of result for next edge  
        if (x != y)  
        {  
        	kv_push(Edge, *result, next_edge);
            Union(subsets, x, y);  
        }  
        // Else discard the next_edge  
    }  
  
    // print the contents of result[] to display the  
    // built MST  
    cout<<"Following are the edges in the constructed MST\n";  
    for (int i = 0; i < result->n; ++i)  
        cout<<result->a[i].vertices[0]<<" -- "<<result->a[i].vertices[1]<<" == "<<result->a[i].weight<<endl;  

    //kv_destroy(*result);
    return;  
}  

int SearchEdgeByVertex(Edges *tree, int v_id)
{
    for (int i = 0; i < tree->n; i++)
    {
        if ( tree->a[i].vertices[0] == v_id || tree->a[i].vertices[1] == v_id)
            return i;
    }
    return -1;
}

pair<int, Edge> get_next_edge(Edges *tree, int cur_edge_id, int vertex_id)
{
    Edge next_edge;
    pair<int,Edge> result;
    int v1,v2;
    for(int i = 0; i < tree->n; i++)
    {
        if( i != cur_edge_id)
        {
            v1 = tree->a[i].vertices[0];
            v2 = tree->a[i].vertices[1];
            if(v1 == vertex_id || v2 == vertex_id )
            {
                int start_pos = (v1 != vertex_id);
                next_edge.vertices[start_pos] = v1;
                next_edge.vertices[!start_pos] = v2;
                next_edge.weight = tree->a[i].weight;
                int result_id = i;
                
                result.first = result_id;
                result.second = next_edge;
                
                return result;
            }
        }

    }

    return  result;
}

//split minimum spanning tree into clusters with size 2-3
void Graph::Split(char *output_fp)
{
    //number of occurence of each vertice 
    int v_occur[V];
    for (int i = 0; i < V; i++)
        v_occur[i] = 0;
    for (int i = 0; i < result->n; i++)
    {
        v_occur[result->a[i].vertices[0]]++;
        v_occur[result->a[i].vertices[1]]++;
    }
    int edge_id = -1, vertex_id = -1;
    //find the id of start verte and pos
    for (int i = 0; i < V; i++)
    {
        int temp_id;
        if ( v_occur[i] == 1 )
        {
            int temp_id = SearchEdgeByVertex(result,i);
            if (edge_id == -1)
            {
                edge_id = temp_id;
                vertex_id = i;
            }
            else if( result->a[temp_id].weight < result->a[edge_id].weight )
            {
                edge_id = temp_id;
                vertex_id = i;
                break;
            }
        }
    }

    int next_edge_id = edge_id;
    //if 0, start from src, else start from dest

    int start_pos = ( result->a[next_edge_id].vertices[0] != vertex_id);

    Edges *sorted_tree = (Edges *)malloc(sizeof(Edges)); // Tnis will store the resultant MST  
    kv_init(*sorted_tree);
    Edge next_edge;
    next_edge.vertices[0] = result->a[next_edge_id].vertices[start_pos];  
    next_edge.vertices[1] = result->a[next_edge_id].vertices[!start_pos];  
    next_edge.weight = result->a[next_edge_id].weight;  
    kv_push(Edge,*sorted_tree,next_edge);

    //sort the tree from start to end
    for(int iter = 0; iter < result->n-1; iter++)
    {
        pair<int,Edge> temp_result = get_next_edge(result, next_edge_id, next_edge.vertices[1]);
        next_edge_id = temp_result.first;
        next_edge = temp_result.second;
        kv_push(Edge,*sorted_tree, next_edge);
    }

   cout<<"Sorted path\n";  
    for (int i = 0; i < sorted_tree->n; ++i)  
        cout<<sorted_tree->a[i].vertices[0]<<" -- "<<sorted_tree->a[i].vertices[1]<<" == "<<sorted_tree->a[i].weight<<endl;  

    cout<<"split path\n";
    ofstream fp;
    fp.open(output_fp);
    int count = 1, start = 0, tree_len = sorted_tree->n;
    for(int i = 0; i < tree_len; i++)
    {
        if (count == 1) {   count++;    }
        else if ( count==2 && i < tree_len-1 && sorted_tree->a[i].weight * threshold < sorted_tree->a[i+1].weight)
        {   count++;    } 
        else
        {
            for(int j = start; j < start+count; j++)
            {
                printf("%d ",sorted_tree->a[i].vertices[0]);
                fp << sorted_tree->a[i].vertices[0] << " ";
            }
            printf("\n");
            fp << "\n";
            start += count;
            count = 0;
        }
    }
    fp.close();
    kv_destroy(*sorted_tree);
    
    return;
}

int minKey(int *key, bool *mstSet, int *degree, int V) 
{ 
    // Initialize min value 
    int min = INT_MAX, min_index; 
 
    for (int v = 0; v < V; v++) 
    {
        if (mstSet[v] == false && key[v] < min && degree[v] < 1) 
        {
            min = key[v];
            min_index = v; 
            degree[v] += 1;
        }
    }
 
    return min_index; 
} 

void PrimMST::RunMST() 
{ 
    //store degree of each edge
    int *degree = (int *)malloc(sizeof(int) * V); 
    // Key values used to pick minimum weight edge in cut 
    int *key = (int *)malloc(sizeof(int) * V);  
    // To represent set of vertices included in MST 
    bool *mstSet = (bool *)malloc(sizeof(bool) * V);
 
    // Initialize all keys as INFINITE 
    for (int i = 0; i < V; i++) 
    {
        key[i] = INT_MAX;
        mstSet[i] = false; 
        degree[i] = 0;
    }
 
    // Always include first 1st vertex in MST. 
    // Make key 0 so that this vertex is picked as first vertex. 
    key[0] = 0; 
    parent[0] = -1; // First node is always root of MST 
 
    int count = 0;
    //V vertices 
    while( count < V - 1 )
    { 

        // Pick the minimum key vertex from the 
        // set of vertices not yet included in MST 
        int u = minKey(key, mstSet,degree,V); 
        // Add the picked vertex to the MST Set 
        mstSet[u] = true; 
 
        // Update key value and parent index of 
        // the adjacent vertices of the picked vertex. 
        // Consider only those vertices which are not 
        // yet included in MST 
        for (int v = 0; v < V; v++) 
 
            // graph[u][v] is non zero only for adjacent vertices of m 
            // mstSet[v] is false for vertices not yet included in MST 
            // Update the key only if graph[u][v] is smaller than key[v] 
            if (dist_matrix[u][v] > 0.0 && mstSet[v] == false && dist_matrix[u][v] < key[v]) 
            {
                parent[v] = u; 
                key[v] = dist_matrix[u][v]; 
            }
    } 
 
    // print the constructed MST 
  //  printMST(parent, dist_matrix);
    printf("result\n") ;
    for (int i = 1; i < V; i++)
        printf("%d-%d\n",parent[i],i);
    delete degree;
    delete key;
    delete mstSet;
} 
