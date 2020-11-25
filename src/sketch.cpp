#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <array>
#include "kvec.h"
#include "creads.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

/*static inline uint64_t hash64_(uint64_t key)
{
	key = (~key + (key << 21)); // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)); // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)); // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31));
	return key;
}
*/
/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param p      minimizers; p->a[i].x is the 2k-bit hash value;
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */

//hash value of minimizer of given reads with quality score larger than threshold

uint32_t sketch_prefix(const char *str,int len,int k)
{
	uint32_t kmer = 0;
	uint64_t mask = (1<<k) - 1;
	for (int i = 0; i < k; i++) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		kmer = (kmer << 2 | c) & mask;           // forward k-mer
		//kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
	}
	return kmer;
}

uint32_t sketch_suffix(const char *str,int len,int k)
{
	uint32_t kmer = 0;
	uint64_t mask = (1<<k) - 1;
	for (int i = len-k; i < len; i++) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		kmer = (kmer << 2 | c) & mask;           // forward k-mer
		//kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
	}
	return kmer;

}


int mm_sketch_single(const char *str, int len,int k)
{
	// printf("\nsingle sequence:");
	// for(int i = 0; i < len; i++)
	// 	printf("%c",str[i]);
	// printf("\n");
	uint64_t mask = (1<<(2*k)) - 1, kmer= 0;//kmer[2] = {0,0};
	int i, l, min=(1<<2*k) - 1;
	//assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		int z;
		kmer = (kmer << 2 | c) & mask;           // forward k-mer
		//kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer

		if (++l >= k && kmer < min) 
		{
			min = kmer;
		}
	}

	return min;
}

int find_minimizer(int *mins, int start, int end, int *pos)
{
	int min = mins[start];
	*pos = start;
//	printf("start:%d end:%d\n",start,end);
	for(int i = start+1; i < end; i++)
	{
//		printf("fm %d,%d\n",i,mins[i]);
		if(mins[i] < min)
		{
			min = mins[i];
			*pos = i;
		} 
	}
	return min;
}

//compute the (w,k) minmizer of the suffix and prefix sbustring,
// s_len <= 2*w+k-1
void mm_sketch_window(const char *str, int len, int w, int k, int s_len, int *result, int *buffer)
{
//	printf("begin sketch window\n");
	//result used as buffer first th
	int l, min=(1<<(2*k)) - 1, kmer = 0, mask = (1<<(2*k))  - 1;
	int result_size = s_len-w-k+2;
	int pre_pos = 0;

//	printf("result size in sketch: %d\n", result_size);
	for( int i = 0, l = 0; i < w+k-1; i++)
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		kmer = (kmer << 2 | c) & mask;		
		if( i >= k-1 )
		{
			buffer[l++] = kmer;
			if ( kmer < min)
			{
				min = kmer;
				pre_pos = i-k+1;
			}
		}

	}
	result[0] = min;
	for (int i = w+k-1, l = 0 ; i < s_len; i++)
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		kmer = (kmer << 2 | c) & mask;
		buffer[i-k+1] = kmer;
		if( pre_pos != l++ && kmer >= min)
		{
			result[l] = min;
		}
		else if (kmer < min)
		{
			min = kmer;
			result[l] = kmer;
			pre_pos = i-k+1;
		}
		else
		{
			result[l] = find_minimizer(buffer, l, l+w, &pre_pos);
			min = result[l];
		}
	}

	int suffix_pos = len-s_len; 
	min=(1<<(2*k)) - 1;
	kmer = 0;
	pre_pos = suffix_pos;

	for( int i = suffix_pos, l = 0; i < suffix_pos+w+k-1; i++)
	{

		int c = seq_nt4_table[(uint8_t)str[i]];
		kmer = (kmer << 2 | c) & mask;		
		if( i >= suffix_pos + k-1 )
		{
			buffer[l++] = kmer;
//			printf("%d kmer: %d\n",l-1, kmer);
			if (kmer < min)
			{
				min = kmer;
				pre_pos = i-k+1;
			}
		}

	}
	result[0+result_size] = min;
//	printf("result %d: %d\n", result_size,result[0+result_size]);
	for (int i = suffix_pos+ w+k-1, l = 1; i < len; i++, l++)
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		kmer = (kmer << 2 | c) & mask;
		buffer[l+w-1] = kmer;
		
	//	printf("kmer %d, pre_pos: %d,%d\n",kmer,pre_pos,suffix_pos+l-1);
		if( pre_pos != suffix_pos+l-1 && kmer >= min)
		{
			result[l+result_size] = min;
		}
		else if (kmer < min)
		{
			result[l+result_size] = kmer;
			min = kmer;
			pre_pos = i-k+1;
		}
		else
		{
			result[l+result_size] = find_minimizer(buffer, l, l+w, &pre_pos);
			pre_pos += suffix_pos;
//			printf("min: %d, pos: %d\n",result[l+result_size],pre_pos);
			min = result[l+result_size];
		}

	}	

	return;
}



/**
void mm_sketch_lh_ori(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, min_pos;
	mm128_t *buf, min = { UINT64_MAX, UINT64_MAX };

	assert(len > 0 && w > 0 && k > 0);
	buf = (mm128_t*)alloca(w * 16);
	memset(buf, 0xff, w * 16);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) 
		{ // not an ambiguous base
			int z;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			if (++l >= k)
				info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		} else l = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k) kv_push(mm128_t, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1) kv_push(mm128_t, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX)
		kv_push(mm128_t, *p, min);
}
**/

int mm_sketch_pos(const char *str, int len,int k)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1<<k) - 1, kmer= 0;//kmer[2] = {0,0};
	int i, l, min=(1<<k) - 1,min_pos = 0;
	assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		int z;
		kmer = (kmer << 2 | c) & mask;           // forward k-mer
		//kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		
		if (++l >= k && kmer < min) 
		{
			min = kmer;
			min_pos = i;
		}
	}

	return min_pos;

}

//return two values, first is hash value of minimizer, second is position
std::array<int, 2> mm_sketch_value_pos(const char *str, int len, int k)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1<<k) - 1, kmer= 0;//kmer[2] = {0,0};
	int i, l, min=(1<<k) - 1,min_pos = 0;
	assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		int z;
		kmer = (kmer << 2 | c) & mask;           // forward k-mer
		//kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		
		if (++l >= k && kmer < min) 
		{
			min = kmer;
			min_pos = i;
		}
	}


	std::array<int, 2> result;
	result[0] = min;
	result[1] = min_pos;
	return result;
}

int mm_sketch_pos2(const char *str, int len,int k)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1<<k) - 1, kmer= 0;//kmer[2] = {0,0};
	int i, l, min=(1<<k) - 1,min_pos = 0;
	assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		int z;
		kmer = (kmer << 2 | c) & mask;           // forward k-mer
		//kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		
		if (++l >= k && kmer < min) 
		{
			min = kmer;
			min_pos = i;
		}
	}
	return min_pos;

}


