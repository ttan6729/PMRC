#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "kvec.h"
#include "creads.h"
#include "cseq.h"
#include "khash.h"
#include <unistd.h>
#include <ios>



using namespace std;

vector<uint32_t> read_numbers; //number of Non-n reads
vector<int> feature_ids;
ofstream info_fp, data_fp;
int **m1; 
float **result_data;
Matrix<int> *minmizers, *buffer;
int **mi_data, **ma_data;
//*mi_norm, *ma_norm;
struct reads_t;
int file_id = 0, num_bucket, seq_len;
int n_nreads; //number of reads contain 'N'
int k,b,w,s_len,result_size;
typedef struct {
	struct reads_t *t;
	long i;
	int tid;
} reads_worker_t; //thread worker

typedef struct reads_t {
	int n_threads;
	int seq_len;
	long n;
	reads_worker_t *w;
	c_seq *seqs;

} reads_t;  //reads for thread work

void mem_usage(double& vm_usage, double& resident_set) {
   vm_usage = 0.0;
   resident_set = 0.0;
   ifstream stat_stream("/proc/self/stat",ios_base::in); //get info from proc
   //create some variables to get info
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;
   unsigned long vsize;
   long rss;
   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
   >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
   >> utime >> stime >> cutime >> cstime >> priority >> nice
   >> O >> itrealvalue >> starttime >> vsize >> rss; 
   stat_stream.close();
   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; //
   vm_usage = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

float Distance(float *v1, float *v2, int len)
{
	float result = 0.0;
	for(int i = 0; i < len; i++)
		result += abs(v1[i] - v2[i]);
	// {
	// 	float value = abs(vectors[id1][i] - vectors[id2][i]);
	// 	if( value < 0)
	// 		value = value*(-1);
	// 	result += value;
	// }
	//	printf("distance %d %d, %.4f\n",id1,id2,result);
	return result;
}


static inline long steal_reads_work(reads_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	return k >= t->n? -1 : k;
}

int element_pos(vector<int> v, int ele)
{
	auto itr = find(v.begin(),v.end(),ele);
	if(itr != v.end())
		return distance(v.begin(),itr);
	return -1;	
}

string strip_slash(const char *c)
{
	string s(c);
	int pos = s.find_last_of("/");
	return s.substr(pos+1);
}


static void process_read(reads_t *t, uint32_t rid, int tid) 
{

	//printf("tid:%d, rid %ld\n",rid);
	int lock=0; //b=14 by default

	c_seq *seq = &t->seqs[rid];
	char *cur_seq = seq->seq;
	// fprintf(stderr, "rid: %ld\n", rid);

	int cntN = 0,cntA = 0, cntT = 0, cntG = 0, cntC = 0;
	for (int i = 0; i < seq_len; ++i) 
	{
		if (cur_seq[i] == 'N')
		{
			__sync_fetch_and_add(&n_nreads,1);
			return;
		}
	}
	// if( rid == 0)
	// 	printf("test paratemeters, seq_len: %d, w: %d, k: %d, s_len: %d, tid %d\n",
	// 		 seq_len, w, k, s_len, tid);
	//printf("tid %d, rid %d\n",tid,rid);

	mm_sketch_window( cur_seq, seq_len, w, k, s_len,minmizers->data[tid],buffer->data[tid] );
	//mm_sketch_single(cur_seq+0, w+k-1,k);
	// printf("string:\n");
	// for(int i = 0; i < seq_len; i++)
	// 	printf("%c,",cur_seq[i]);
	// printf("\nassert:\n");
	for(int i  = 0; i < result_size; i++)
	{
		__sync_fetch_and_add(&mi_data[file_id][minmizers->data[tid][i]],1);
	 	//assert(minmizers->data[tid][i]==mm_sketch_single(cur_seq+i, w+k-1,k));
	 	//printf("result: %d exptected: %d\n",minmizers->data[tid][i], mm_sketch_single(cur_seq+i, w+k-1,k));
	}

	// printf("suffix string:\n");
	// for(int i = seq_len-s_len; i < seq_len; i++)
	// 	printf("%c,",cur_seq[i]);
	// printf("\n");
	for(int i  = result_size, l = 0; i < 2*result_size; i++,l++)
	{
		__sync_fetch_and_add(&mi_data[file_id][num_bucket+minmizers->data[tid][i]],1);	
		//assert(minmizers->data[tid][i]==mm_sketch_single(cur_seq+seq_len-s_len+l, w+k-1,k));
	 	//printf("result %d: %d exptected: %d\n",i,minmizers->data[tid][i],mm_sketch_single(cur_seq+seq_len-s_len+l, w+k-1,k));
	}



//	lock = mm_sketch_single(cur_seq, seq_len,k);//mm_sketch_single(cur_seq, seq_len,b);
	//__sync_fetch_and_add(&m1[file_id][lock],1);
	return;
}

static void *reads_worker(void *data) 
{
	// int test_number =1000;
	// printf("test number:%d\n",test_number);

	reads_worker_t *w = (reads_worker_t*)data;
	reads_t *t = w->t;
	long i;
	for (;;) 
	{
		i = __sync_fetch_and_add(&w->i, t->n_threads);
	//	if (i > test_number) return NULL;
		if (i >= t->n) break;
		process_read(t, i, w->tid);
	}
	while ((i = steal_reads_work(t)) >= 0)
		process_read(t, i,w->tid);
	pthread_exit(0);
}

void process_file(const char *file_name, int length,reads_t t)
{

	int i;
	seq_file *fp = seq_open(file_name);
	//if (fp == 0) return; 
	int n_seq = 0;
	c_seq *seqs = seq_init(fp,&n_seq,&length);

//	info_fp << strip_slash(file_name) << " " << n_seq << "\n";
	t.seqs = seqs;
	t.n = n_seq;
	n_nreads = 0;
//	printf("read length:%d, seq number:%d,thread number:%d\n",length,n_seq,t.n_threads);
	pthread_t *tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) 
		t.w[i].t = &t, t.w[i].i = i, t.w[i].tid = i; //
	// for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_worker, &t.w[i]);

	for (i = 0; i < n_threads; ++i)
		pthread_create(&tid[i], 0, reads_worker, &t.w[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], NULL);

	read_numbers.push_back(n_seq - n_nreads);
	for(int i = 0; i < n_seq; i++)
		free(seqs[i].seq);
	free(seqs);
	seq_close(fp);

	return;	
}


void init_file()
{

}
template <typename T>
void print_data( T **data, int n, int m)
{
	for( int i = 0; i < n; i++ )
	{
		for( int j = 0; j < m; j++)
			cout << data[i][j] << " ";
		cout << endl;
	}

}
float **pre_process(vector<string> file_list, char *dir, int _k, int selected_number, char *matrix_fp)
{
	printf("begin process, thread number %d\n",n_threads);
	k = _k;
	w = k;
	b = 2 * k;
	reads_t t;
	pthread_t *tid;

	char info_name[1024];
	sprintf(info_name,"%s/info",dir);
	info_fp.open(info_name,std::fstream::in | std::fstream::out | std::fstream::app);
	for(int i = 0; i < file_list.size(); i++)
		info_fp << strip_slash(file_list[i].c_str()) << "\n";
	info_fp.close();
	
	t.n_threads = n_threads;//n_threads;
	t.w = (reads_worker_t*)alloca(n_threads * sizeof(reads_worker_t));

    num_bucket = 1<<b;
	
	mi_data = (int**)malloc(sizeof(int *) * file_list.size());
	for(int i = 0; i < file_list.size(); i++)
	{
		mi_data[i] = (int *)malloc(sizeof(int)*num_bucket*2);
		memset(mi_data[i],0,sizeof(mi_data[i][0])/sizeof(char)*num_bucket*2);
	}

	seq_len = 0;
	for(int i = 0; i < file_list.size(); i++)
	{
		ifstream file;
	   	file.open(file_list[i],ios::in); 
	   	if (file.is_open())
	   	{   //checking whether the file is open
	      string line;
	      getline(file, line);
	      getline(file, line);
	      seq_len = strlen(line.c_str());

	      if (i == 0)	
	      {
	      	s_len = w+k-1+w-1;
		  	if (s_len >= seq_len/2)
		  		s_len = seq_len/2;
	      	t.seq_len = seq_len;
	      	result_size = s_len-w-k+2;
			minmizers = new Matrix<int>(n_threads,2*result_size);
			buffer = new Matrix<int>(n_threads,s_len-k+1);
			printf("test parameter:k:%d w:%d s_len:%d\n",k,w,s_len);
	      }
	      else if ( t.seq_len != seq_len)
	      {
	      	printf("error, the length of dataset doesn't equal\n");
	      	exit(0);
	      }
	      //file.close(); //close the file object.
	      process_file(file_list[i].c_str(),seq_len,t);
	     // printf("test for first file\n");
	      file_id++;
	    }
	    else 
	    {
	    	printf("failed to open %s\n",file_list[i].c_str());
	    } 
	}
	//return NULL;

	//normalization
	//mi_norm = new Matrix<double>(file_list.size(),num_bucket);
	//ma_norm = new Matrix<double>(file_list.size(),num_bucket);
	Matrix<double> *norm_data = new Matrix<double>(file_list.size(),2*num_bucket);
	for(int i = 0; i < file_list.size(); i++)
	{
		for(int j = 0; j < 2*num_bucket; j++)
		{
			norm_data->data[i][j] = 100.0*mi_data[i][j]/read_numbers[i];
		//	norm_data[i+num_bucket][j] = 100.0*ma_data[i][j]/read_numbers[i];		
			// mi_norm->data[i][j] = 100.0*mi_data[i][j]/read_numbers[i];
			// ma_norm->data[i][j] = 100.0*ma_data[i][j]/read_numbers[i];
		}
	}

	vector<double> vars;
	for(int i = 0; i < 2*num_bucket; i++)
		vars.push_back(norm_data->col_var(i));
	vector<double> vars2(vars);
	sort(vars2.begin(), vars2.end());

	result_data = (float**)malloc(sizeof(float *) * file_list.size());
	for(int i = 0; i < file_list.size(); i++)
		result_data[i] = (float *)malloc(sizeof(float) *selected_number);
	ofstream vector_fp;
	int mi_count = 0;
	for(int i = 0; i < selected_number; i++)
	{
		double value = vars2[vars2.size()-i-1];
		auto itr = find(vars.begin(),vars.end(),value);
		int index = distance(vars.begin(),itr);
		feature_ids.push_back(index);
		if ( index < num_bucket)
			mi_count++;
		//printf("feature id:%d,vars %.2f,vars check %.2f\n",index,value,norm_data->col_var(index));
		for(int j = 0; j < file_list.size(); j++)
			result_data[j][i] = norm_data->data[j][index];

	}
	printf("%d features from minimizer, %d from maximizer\n",mi_count,selected_number-mi_count);
	free(mi_data);

	Matrix<float> *dist_matrix =  new Matrix<float>(file_list.size(),file_list.size());
	for(int i = 0; i < file_list.size(); i++)
		data_fp << "," << strip_slash(file_list[i].c_str()) ;
	data_fp << "\n";
	for(int i = 0; i < file_list.size()-1; i++)
	{
		for(int j = i+1; j < file_list.size(); j++)
		{
			dist_matrix->data[i][j] = Distance(result_data[i],result_data[j],selected_number); 
			dist_matrix->data[j][i] = dist_matrix->data[i][j];
		}
	}
	for(int i = 0; i < file_list.size(); i++)
		dist_matrix->data[i][i] = 0.0;
	
	if(matrix_fp != NULL)
	{
		printf("write distance matrix to %s\n",matrix_fp);
		data_fp.open(matrix_fp);
		for(int i = 0; i < file_list.size(); i++)
			data_fp << "," << strip_slash(file_list[i].c_str());
		data_fp << "\n";

		for(int i = 0; i < file_list.size(); i++)
		{
			dist_matrix->data[i][i] = 0.0;
			data_fp << strip_slash(file_list[i].c_str());
			for(int j = 0; j < file_list.size(); j++)
				data_fp << "," << dist_matrix->data[i][j];
			data_fp << "\n";
		}
		data_fp.close();
	}


	return dist_matrix->data;
}