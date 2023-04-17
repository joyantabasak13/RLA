// Implement ukkonen's aproximate string mapping algorithm
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iostream>
#include <map>
#include <pthread.h>
#include <set>
#include <string>
#include <sys/time.h>
#include <random>
#include <tuple>
#include <time.h>
#include <unistd.h>
#include <utility>
#include <vector>
#include <mutex>
#include <queue>

// #include <bits/stdc++.h>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

using namespace std;

int threshold = 99;
int matArr[50][50] = {0};

// helps edit distance calculation in calculateBasicED()
int calculateBasicED2(string &str1, string &str2, int threshRem)
{
	int row, col, i, j;
	row = str1.length() + 1;
	col = str2.length() + 1;

	// matArr.resize(row);
	// for (i = 0; i < row; ++i)
	// 	matArr[i].resize(col, 0);

	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
		{
			if (i == 0)
				matArr[i][j] = j;
			else if (j == 0)
				matArr[i][j] = i;
			else
			{
				if ((int)str1[i - 1] == (int)str2[j - 1])
					matArr[i][j] = min(min(matArr[i - 1][j - 1], matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
				else
					matArr[i][j] = min(min(matArr[i - 1][j - 1] + 1, matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
			}

			if ((row - col) == (i - j) && (matArr[i][j] > threshRem))
			{
				return threshold + 1;
			}
		}
	}
	
	return (matArr[row - 1][col - 1]);
}

// calculates edit distance between two string
// takes two strings and a threshold value as input
// returns global variable threshold + 1 if distance exceeds theshold param
// returns edit distance
// core mechanism is a DP algo
int calculateBasicED(string &str1, string &str2, int threshRem)
{
	int dist = threshRem;

	if (abs((int)(str1.length() - str2.length())) > dist)
		return threshold + 1;
	else if (str1.compare(str2) == 0)
		return 0;
	else if (dist == 0)
		return threshold + 1;
	else if ((2 * dist + 1) >= max(str1.length(), str2.length()))
		return calculateBasicED2(str1, str2, dist);
	else
	{
		
		string s1, s2;
		int row, col, diagonal;
		int i, j;

		if (str1.length() > str2.length())
		{
			s1 = str2;
			s2 = str1;
		}
		else
		{
			s1 = str1;
			s2 = str2;
		}

		row = s1.length() + 1;
		col = 2 * dist + 1;
		diagonal = dist + s2.length() - s1.length();

		// matArr.resize(row);
		// for (i = 0; i < row; ++i)
		// 	matArr[i].resize(col, 0);

		// if(procID == 1 && checkTemp == 3164)
		//	cout << str1 << " -- " << str2 << " rt " << dist << endl;

		for (i = 0; i < dist + 1; i++)
		{
			for (j = dist - i; j < col; j++)
			{
				if (i == 0)
					matArr[i][j] = j - dist;
				else if (j == (dist - i))
					matArr[i][j] = matArr[i - 1][j + 1] + 1;
				else if (j != (col - 1))
				{
					if ((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
						matArr[i][j] = min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] = min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				else
				{
					if ((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
						matArr[i][j] = min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}

				if ((j == diagonal) && matArr[i][j] > dist)
					return threshold + 1;
			}
		}

		for (i = dist + 1; i < s2.length() - dist + 1; i++)
		{
			for (j = 0; j < col; j++)
			{
				if (j == 0)
				{
					if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] = min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else if (j != (col - 1))
				{
					if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] = min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] = min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				else
				{
					if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] = min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}
				if ((j == diagonal) && (matArr[i][j] > dist))
					return threshold + 1;
			}
		}

		for (i = s2.length() - dist + 1; i < row; i++)
		{
			for (j = 0; j < col - i + s2.length() - dist; j++)
			{
				if (j == 0)
				{
					if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] = min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else
				{
					if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] = min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] = min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				if ((j == diagonal) && (matArr[i][j] > dist))
					return threshold + 1;
			}
		}

		// if(procID == 1 && checkTemp == 3164)
		// cout << str1 << " -- " << str2 << " hukjhk " << matArr[row - 1][diagonal] << endl;

		return matArr[row - 1][diagonal];
	}
}

ptrdiff_t algo_8(const std::string &A, ptrdiff_t m, const std::string &B, ptrdiff_t n, ptrdiff_t inf, ptrdiff_t k, ptrdiff_t p, size_t num_p, ptrdiff_t *f)
{
	ptrdiff_t t = f[num_p * (k + m) + (p)] + 1; // f(k,p-1) + 1
	ptrdiff_t t1 = -inf;
	if (k + m >= 1 && p >= 0)
	{
		t1 = f[num_p * (k + m - 1) + (p)];
	}

	ptrdiff_t t2 = -inf;
	if (k < n && p < inf + 2)
	{
		t2 = f[num_p * (k + m + 1) + (p)] + 1; // f(k+1,p-1) + 1
	}

	// check whether a(t)a(t+1) = b(k+t+1)b(k+t), but *only* if those
	// ranges are valid. IOW, if:
	//
	//   - t > 0
	//   - t + 1 <= m
	//   - k + t > 0
	//   - k + t + 1 <= n
	ptrdiff_t t3 = -inf;
	if (t > 0 &&
		t + 1 <= m &&
		k + t > 0 &&
		k + t + 1 <= n)
	{
		if (A[t - 1] == B[k + t] && A[t] == B[k + t - 1])
		{
			t3 = t + 1;
		}
	}
	// t := max(t,t1,t2,t3)
	if (t < t1)
	{
		t = t1;
	}
	if (t < t2)
	{
		t = t2;
	}
	if (t < t3)
	{
		t = t3;
	}
	// while a(t+1) = b(t+1+k) do t := t + 1
	while (t + 1 <= m &&
		   t + k + 1 <= n &&
		   A[t] == B[t + k])
		t++;

	if (t > m || t + k > n)
	{
		t = inf;
	}
	return t;
}

ptrdiff_t ukkonen(const std::string &A, const std::string &B, std::size_t D)
{
	ptrdiff_t m = A.length();
	ptrdiff_t n = B.length();
	ptrdiff_t inf = max(m, n); // |A,B| <= inf

	// Allocating too much space, here. This implementation won't satisfy the
	// space bounds.

	// To index into `f' in terms of (i,j), -m <= i <= n,
	// -1 <=j <= inf, do f[i+m][j+1]
	ptrdiff_t f[(m + n + 1) * (inf + 2)];
	for (ptrdiff_t i = 0; i < m + n + 1; ++i)
	{
		for (ptrdiff_t j = 0; j < inf + 2; ++j)
		{
			f[i * (inf + 2) + j] = -inf - 1;
		}
	}

	// Initialize f: f(k,|k|-1) = |k|-1, if k < 0...
	for (ptrdiff_t k = -1; k >= -m; --k)
	{
		f[(k + m) * (inf + 2) - k] = -k - 1;
	}
	for (ptrdiff_t k = 0; k <= n; ++k)
	{
		f[(k + m) * (inf + 2) + k] = -1;
	}

	ptrdiff_t p = -1;
	ptrdiff_t r = p - min(m, n);
	while (f[n * (inf + 2) + p + 1] != m)
	{
		p = p + 1;
		r = r + 1;
		if (r <= 0)
		{
			for (ptrdiff_t k = -p; k <= p; ++k)
			{
				// f(k,p)
				f[(inf + 2) * (k + m) + p + 1] = algo_8(A, m, B, n, inf, k, p, inf + 2, (ptrdiff_t *)f);
			}
		}
		else
		{
			for (ptrdiff_t k = max(-m, -p); k <= -r; ++k)
			{
				// f(k,p)
				f[(inf + 2) * (k + m) * p + 1] = algo_8(A, m, B, n, inf, k, p, inf + 2, (ptrdiff_t *)f);
			}
			for (ptrdiff_t k = r; k <= max(n, p); ++k)
			{
				// f(k,p)
				f[(inf + 2) * (k + m) + p + 1] = algo_8(A, m, B, n, inf, k, p, inf + 2, (ptrdiff_t *)f);
			}
		}
		if(p>D) {
			return threshold+1;
		}
	}
	return p;
}

double jaro_distance(string& s1, string& s2)
{
	// Length of two strings
	int len1 = s1.length();
	int len2 = s2.length();

	if (len1 == 0 || len2 == 0)
		return 0.0;

	// Maximum distance upto which matching is allowed
	int max_dist = floor(max(len1, len2) / 2) - 1;

	// Count of matches
	int match = 0;

	// Hash for matches
	std::vector<int> hash_s1(s1.length(),0);
	std::vector<int> hash_s2(s2.length(),0);

	// Traverse through the first string
	for (int i = 0; i < len1; i++) {

		// Check if there is any matches
		for (int j = max(0, i - max_dist);
			j < min(len2, i + max_dist + 1); j++)
			// If there is a match
			if (s1[i] == s2[j] && hash_s2[j] == 0) {
				hash_s1[i] = 1;
				hash_s2[j] = 1;
				match++;
				break;
			}
	}

	// If there is no match
	if (match == 0)
		return 0.0;

	// Number of transpositions
	double t = 0;

	int point = 0;

	// Count number of occurrences
	// where two characters match but
	// there is a third matched character
	// in between the indices
	for (int i = 0; i < len1; i++)
		if (hash_s1[i]) {

			// Find the next matched character
			// in second string
			while (hash_s2[point] == 0)
				point++;

			if (s1[i] != s2[point++])
				t++;
		}

	t /= 2;

	// Return the Jaro Similarity
	return (((double)match) / ((double)len1)
			+ ((double)match) / ((double)len2)
			+ ((double)match - t) / ((double)match))
		/ 3.0;
}




// double calculateBasicJaroWinklerDist(string& str1,string& str2,double threshRem) {

// 	//cout << "Strings\t" << str1 << "," << str2 << "\n";
// 	int min_string_len = (str1.size() - str2.size());
	
// 	if (min_string_len < 0){
// 		min_string_len = -min_string_len;
// 	}
// 	if (min_string_len > 1) {
// 		return 0.20;
// 	}

// 	double jaro_dist = jaro_distance(str1, str2);

// 	//cout << "Jaro Distance" << jaro_dist << "\n";
// 	// If the jaro Similarity is above a threshold
// 	if (jaro_dist > threshRem) {
// 	    //	cout << "Inside here" << "\n";
// 		// Find the length of common prefix
// 		int prefix = 0;

// 		for (int i = 0;
// 			i < min(str1.length(), str2.length()); i++) {
// 			// If the characters match
// 			if (str1[i] == str2[i])
// 				prefix++;

// 			// Else break
// 			else
// 				break;
// 		}

// 		// Maximum of 4 characters are allowed in prefix
// 		prefix = min(4, prefix);

// 		// Calculate jaro winkler Similarity
// 		jaro_dist += 0.1 * prefix * (1 - jaro_dist);
// 	}
// 	//cout << "Jaro_winkler" <<  jaro_dist << "\n";
// 	return jaro_dist;

// }

// double calculateJaroWinklerDist(vector<string>& a, vector<string>& b, vector<vector<int> >& compareAtt,int threshRem) {
	
// 	int ind, a_ind, b_ind;
// 	double temp, w;
// 	string s1, s2;

// 	//cout << "Inside Hausdorff";

// 	w 	= 0;

// 	int al_ind = atoi(a[a.size() - 1].c_str()), bl_ind = atoi(b[b.size() - 1].c_str());
// 	for (int g = 0; g < compareAtt.size(); g++)
// 	{
// 		ind 	= compareAtt[g][0];
// 	//	cout << ind << "ind";
// 		a_ind	= indexDatasetArr[al_ind][ind];
// 		b_ind	= indexDatasetArr[bl_ind][ind];

// 		if (a_ind == -1 || b_ind == -1 || a[a_ind].length() == 0 || b[b_ind].length() == 0)
// 			continue;

// 		c2++;
// 		s1 		= a[a_ind];
// 		s2 		= b[b_ind];

		
// 		//  temp 	= calculatePosHausdorffDistance(s1, s2,1) * weightArr[ind];
// 		//cout << "S1,S2" << s1 << "," << s2 << "\n";
// 		temp 	= calculateBasicJaroWinklerDist(s1, s2,threshold) * weightArr[ind];
// 		//  temp = calculateBasicHausdorffDistance(s1,s2,0.55);
// 		w 			+= temp;
		
// 		//cout << "W" << w << "\n";
// 		//h_distance.insert({string_pair,w});
		
		
// 		if (((w + 1) / 2) < threshold )
// 			break;
		
// 	}

// 	return (w / compareAtt.size());
// }


int main(int argc, char **argv)
{
	string A = "Joyanta";
	string B = "Joyanta";

	int dist = calculateBasicED(A, B, 5);
	double jaro_dist = jaro_distance(A,B);
	cout << "Edit distance: " << dist << endl;
	cout<< "Jaro Similarity: " << jaro_dist << endl;

	return 0;
}



