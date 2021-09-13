// testing1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#define MAX_LEN 256
#define SAMPLE_SIZE 320
#define past_samples 12
#define PI 3.14
char vowels[5] = {'a','e','i','o','u'};
double sin_w[12] = {2.552145, 3.998621, 5.240951, 6.194559, 6.794523, 6.999998, 6.796995, 6.199335, 5.247705, 4.006894, 2.561374, 1.009556};


/* this method will do the verification of capstral coefficients (Ci) on the test data */
bool coefficients_verfication(){
	FILE* file = fopen("resources/test.txt", "r");
	double arr[SAMPLE_SIZE], a[SAMPLE_SIZE+1], r[SAMPLE_SIZE+1], c[SAMPLE_SIZE+1], temp[SAMPLE_SIZE+1][SAMPLE_SIZE+1], temp_var, r_i;
	int i=0,j, k, m;
	bool is_verf_true = true;
	// -1 to allow room for NULL terminator for really long string
    while (fscanf(file,"%lf\n",&temp_var)!= EOF)
    {	
		arr[i++] = temp_var;
		//printf("%ld \n", temp_var);
    }
	fclose(file);
	
	//-------checking Ri----------------------
	//printf("-------checking Ri----------------------\n");
	for(k=0; k<past_samples+1; k++)		//since we need values from r[0] to r[past_samples]
	{
		r_i = 0;
		for(m=0; m<=(SAMPLE_SIZE-1-k); m++)
			r_i += (arr[m]*arr[m+k]);

		r[k] = r_i;
		//printf("%lf\n", r[k]);
	}

	//-------checking Ai----------------------
	//printf("-------checking Ai----------------------\n");
	a[0] = r[0];
	if(a[0]==0)
	{
		return false;
	}

	double sum, temp2;
	for(i=1; i<=past_samples; i++)
	{
		sum = 0;
		for(j=1; j<=(i-1); j++)
			sum += temp[i-1][j]*r[i-j];

		temp2 = (r[i]-sum)/a[i-1];
		
		temp[i][i] = temp2;

		for(j=1; j<=i-1; j++)
			temp[i][j] = temp[i-1][j] - (temp2*temp[i-1][i-j]);

		a[i] = (1 - temp2*temp2) * a[i-1];
	}

	for(i=1; i<=past_samples; i++)
	{
		a[i] = temp[past_samples][i];
		//printf("%lf\n", a[i]);
	}

	//-------checking Ci----------------------
	double temp3;
	//printf("-------checking Ci----------------------\n");
	for(i=1; i<=past_samples; i++)
	{
		sum = 0;
		for(k=1; k<=(i-1); k++)
			{
				temp3 = ((double)k/(double)i)*c[k]*a[i-k];
				//printf("\nnum=%Lf",num);
				sum += temp3;
			}

		c[i] = a[i]+sum;	
		//printf("%lf\n", c[i]);
	}

	return true;
}

/*this method will calculate the ste for each of the sample mentioned by i */
double find_energy(double *a, int i){
	double ste = 0;
	for(int t = i; t < i+100; t++){
		ste = ste + a[t]*a[t];
	}
	ste = ste/100;
	return ste;
}

/* this method will help us in finding the start marker for each of the file*/
int find_start_marker(char vowel, int fileno){
	char src[128];
	_snprintf(src, sizeof(char) * 128, "%s%c%s%d%s", "resources/processed/normalization/", vowel, "_", fileno, ".txt");
	double abs_max = -1000000, temp_var, a[100000];
	int i,k, start_index=0,n=0;
	char buffer[MAX_LEN];

	FILE* file = fopen (src, "r");
	while (fscanf(file,"%lf\n",&temp_var)!= EOF)
    {	
		a[n++] = temp_var;
    }
    fclose(file);

	//calculating abs max
	/*
	for(i=0; i<n; i++){

		if(fabs(a[i])>abs_max){
			abs_max =fabs(a[i]);
		} 
	}

	for(i=0; i<n; i++){
		a[i] = ((a[i]*5000)/ abs_max);
	}
	*/

	double ts = 0;
	for(int i=0; i<5; i++){
		ts = ts + find_energy(a,i*100);
	}
	ts = ts /5.0;

	double THRESHOLD_SILENCE = 10* ts;
	//printf("Threshold: %lf \n", ts);

	for(i=1000; i<n-1000; i++){
		if(find_energy(a,i)> THRESHOLD_SILENCE){
			int count=0;
			for(k=1; k<=10; k++){
				if(find_energy(a,i+(k*100))>THRESHOLD_SILENCE){
					count++;
				}
			}
			if(count>=8){
				return i;
			}
		}
	}
	return -1;
}

/* this method will calculate ri, ai and then ci. and finally store ci in the file. */
bool calc_capstral_coeff(char vowel, int fileno, int strt){
	char src[128], dest[128];
	_snprintf(src, sizeof(char) * 128, "%s%c%s%d%s", "resources/processed/normalization/", vowel, "_", fileno, ".txt");
	_snprintf(dest, sizeof(char) * 128, "%s%c%s%d%s", "resources/processed/cc_calc/", vowel, "_", fileno, ".txt");
	
	double data[320*5], temp_var;
	int cntr=0;

	//reading the file from start marker in an array---------------------------------------------------
	FILE* file = fopen(src, "r");
	while (fscanf(file,"%lf\n",&temp_var)!= EOF)
    {	
		if(cntr< strt){
			cntr++;
			continue;
		}
		if(cntr== strt+320*5){
		    break;
		}
		data[cntr-strt] = temp_var;
		cntr++;
    }
    fclose(file);

	//writing the capstral coeff to the file---------------------------------------------------
	file = fopen (dest, "w");
	for(int f=0; f<5; f++){ //for 5 frames

		double arr[SAMPLE_SIZE], a[SAMPLE_SIZE+1], r[SAMPLE_SIZE+1], c[SAMPLE_SIZE+1], temp[SAMPLE_SIZE+1][SAMPLE_SIZE+1], temp_var, r_i;
		int i,j,k,m;
		for(int s=0; s< SAMPLE_SIZE; s++){
			arr[s] = data[(f*SAMPLE_SIZE) + s];
		}
		//-------checking Ri----------------------
		//printf("-------checking Ri----------------------\n");
		for(k=0; k<past_samples+1; k++)		//since we need values from r[0] to r[past_samples]
		{
			r_i = 0;
			for(m=0; m<=(SAMPLE_SIZE-1-k); m++)
				r_i += (arr[m]*arr[m+k]);

			r[k] = r_i;
			//printf("%lf\n", r[k]);
		}

		//-------checking Ai----------------------
		//printf("-------checking Ai----------------------\n");
		a[0] = r[0];
		if(a[0]==0)
		{
			return false;
		}

		double sum, temp2;
		for(i=1; i<=past_samples; i++)
		{
			sum = 0;
			for(j=1; j<=(i-1); j++)
				sum += temp[i-1][j]*r[i-j];

			temp2 = (r[i]-sum)/a[i-1];
			temp[i][i] = temp2;

			for(j=1; j<=i-1; j++)
				temp[i][j] = temp[i-1][j] - (temp2*temp[i-1][i-j]);

			a[i] = (1 - temp2*temp2) * a[i-1];
		}

		for(i=1; i<=past_samples; i++)
		{
			a[i] = temp[past_samples][i];
			//printf("%lf\n", a[i]);
		}

		//-------checking Ci----------------------
		double temp3;

		//printf("-------checking Ci----------------------\n");
		double w = 0,z;
		for(i=1; i<=past_samples; i++)
		{
			sum = 0;
			for(k=1; k<=(i-1); k++)
				{
					temp3 = ((double)k/(double)i)*c[k]*a[i-k];
					//printf("\nnum=%Lf",num);
					sum += temp3;
				}
			w = 1.0 + 6.0 * sin((PI * i)/12.0);
			c[i] = (a[i]+sum);
			//z=c[i];
			//c[i] = c[i] * w;
			//printf("%lf- %lf- %lf \n",w,z, c[i]);
			//printf("%lf \n",c[i]);
			fprintf(file, "%lf\n", c[i]);
		}
	}
	fclose(file);

	return true;
} 


/* main fuction which is used for dc shift and normlalization */
int DC_Shift_and_Normalization(char vowel, int fileno){
	char src[128], dest[128];
	_snprintf(src, sizeof(char) * 128, "%s%c%s%d%s", "resources/recordings/214101048_", vowel, "_", fileno, ".txt");
	_snprintf(dest, sizeof(char) * 128, "%s%c%s%d%s", "resources/processed/normalization/", vowel, "_", fileno, ".txt");

	char buffer[MAX_LEN];
	FILE* file = fopen (src, "r");

	double abs_max = -1000000, dc_shift = 0, x=0;
	int count=0;

	int i=0;
	while (fgets(buffer, MAX_LEN - 1, file))
    {	
		i++;
		//ignoring first 10 frame as it may contain noise.
		if(i<1000){
			continue;
		}
        // Remove trailing newline
        buffer[strcspn(buffer, "\n")] = 0;
		
		x = (double) atoi(buffer);
		
		//calcuating abs max  which will help in normalization of the data.
		if(fabs(x)>abs_max){
			abs_max = fabs(x);
		}    

		//for dc shift
		count++;
		dc_shift += x;
    }

	//dc_shift = dc_shift/ count;
	dc_shift = 0;
    fclose(file);

	//printf("DC Shift is %lf and max value of array is %lf for normalization. \n",dc_shift,abs_max);
	FILE* fp1 = fopen (src, "r");
	FILE* fp2 = fopen (dest, "w");

	i=0;x=0;
	while (fgets(buffer, MAX_LEN - 1, fp1))
	{
		i++;
		
		//ignoring first 10 frame as it may contain noise.
		if(i<1000){
			continue;
		}

        // Remove trailing newline
        buffer[strcspn(buffer, "\n")] = 0;
		
		x = (double) atoi(buffer);

		x = (((x - dc_shift)*5000)/ abs_max);

		fprintf(fp2, "%lf\n", x);
	
	}

	fclose(fp1);
	fclose(fp2);
	
	return 1;
}

/*helper function of predict_vowel class to calculate tokura distance..*/
double tokura_distance(double c_t[],double c_r[], int frame)
{
	double w[12] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
	double distance=0,x=0;
	int i=0;
	for(i=(frame*past_samples);i<((frame+1)*past_samples);i++)
	{
		x=((c_t[i]*sin_w[i-(frame*past_samples)])-(c_r[i]*sin_w[i-(frame*past_samples)]));
		distance+=(w[i-(frame*past_samples)]*x*x);
	}
	//printf("%lf ", distance);
	return distance;
}

/* main function for predicting the vowel class */
char predict_vowel(char vowel,int fileno){
	double c_t[5*past_samples],c_r[5*past_samples];
	double dist[5], temp_var;
	char test[128], reference[128];
	int cntr;
	FILE *tf_p, *rf_p;
	_snprintf(test, sizeof(char) * 128, "%s%c%s%d%s", "resources/processed/cc_calc/", vowel, "_", fileno, ".txt");
	tf_p = fopen(test, "r");
	cntr=0;
	while (fscanf(tf_p,"%lf\n",&temp_var)!= EOF)
	{	

		c_t[cntr++] = temp_var;
	}
	fclose(tf_p);

	for(int rf=0; rf<5; rf++){ // checking with all the 5 vowels
		_snprintf(reference, sizeof(char) * 128, "%s%c%s%s", "resources/processed/training/", vowels[rf], "_train",".txt");
		rf_p = fopen(reference, "r");
		cntr=0;
		while (fscanf(rf_p,"%lf\n",&temp_var)!= EOF)
		{	
			c_r[cntr++] = temp_var;
		}
		fclose(rf_p);
		double avg_dist=0;
		for(int f=0; f<5; f++){ //for each frame
			avg_dist = avg_dist + tokura_distance(c_t,c_r,f);
		}
		avg_dist = avg_dist/5.0;
		dist[rf] = avg_dist;
	}

	//predicting vowel----
	double minm = dist[4];
	int index = 4;
	for(int v=0; v<5; v++){
		if(dist[v]< minm){
			minm = dist[v];
			index = v;
		}
	}
	
	printf("\n\nTokura distances: ");
	for(int i=0; i<5; i++){
		printf("%lf ", dist[i]);
	}
	printf("\n");
	return vowels[index];
}

int _tmain(int argc, _TCHAR* argv[])
{
	char buffer[MAX_LEN];
	int n=0, x,i=0;


	//-----------step 1: Verification of capstral coefficient
	
	bool valid = coefficients_verfication();

	if(valid){
		printf( "Ai Ri Ci coefficients verification passed. Great !! \n" );
	}else{
		printf( "Ai Ri Ci coefficients verification failed. Try again !! \n" );
	}

	printf( "Loading the model.... and starting the prediction!! \n" );

	/*------------ step 2: Get features of Reference recordings. */
	char fname[128];
	
	for(int i=0; i<5; i++){
		for(int j=1; j<=20; j++){
			char vowel = vowels[i];
			int fileno = j;
			DC_Shift_and_Normalization(vowel,fileno);
		}
	}
	
	//---------------------finding start marker
	int strt_mar[5][20];
	for(int i=0; i<5; i++){
		for(int j=1; j<=20; j++){
			char vowel = vowels[i];
			int fileno = j;
			int start_marker = find_start_marker(vowel,fileno);
			strt_mar[i][j-1] = start_marker;
			//printf("%c-%d-%d\n", vowel, fileno, start_marker);
		}
	}

	/*-----------------step 3: Get features of training & test recordings. */
	for(int i=0; i<5; i++){
		for(int j=1; j<=20; j++){
			char vowel = vowels[i];
			int fileno = j;
			calc_capstral_coeff(vowel, fileno, strt_mar[i][j-1]);
			//printf("%c-%d\n", vowel, fileno);
		}
	}

	// Training features ------------------------------------------------------------------
	for(int i=0; i<5; i++){
		char vowel = vowels[i];
		char dest[128];
		double featres[5*past_samples], temp_var;
		for(int k=0; k< 5*past_samples; k++){
			featres[k] = 0.0;
		}
		_snprintf(dest, sizeof(char) * 128, "%s%c%s%s", "resources/processed/training/", vowel, "_train",".txt");
		FILE* fp1 = fopen(dest, "w");

		for(int j=1; j<=10; j++){
			int fileno = j;
			char src[128];
			_snprintf(src, sizeof(char) * 128, "%s%c%s%d%s", "resources/processed/cc_calc/", vowel, "_", fileno, ".txt");
			//printf("vowel-%c- File-%d\n", vowel,fileno);
			int cntr=0;
			FILE* fp2 = fopen(src, "r");
			while (fscanf(fp2,"%lf\n",&temp_var)!= EOF)
			{	
				featres[cntr] = featres[cntr] + temp_var;
				cntr++;
			}
			fclose(fp2);
			
		}

		for(int t=0; t<5*past_samples; t++){
			featres[t] = featres[t]/10.0;
			fprintf(fp1, "%lf\n", featres[t]);
		}
		fclose(fp1);
	}

	/* step 4: Find similarity and do the predictions.Total 50 test files: 10 for each vowels.*/ 
	double count_corr_pred = 0;
	for(int v=0; v<5; v++){  //vowel no.
		for(int f=11; f<=20; f++){ //test file no.
			char vowel = vowels[v];
			int fileno = f;
			char pred = predict_vowel(vowel,fileno);

			if(pred == vowel){
				count_corr_pred++;
			}
			printf("Actual Class:[%c , %d] \t Predicted Class:%c \n", vowel,fileno, pred);
		}
	}

	printf("\n\nFinal accuracy is %lf Percent \n", (count_corr_pred/50.0)*100.0);


	getchar();
	return 0;
}