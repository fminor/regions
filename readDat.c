/************************************************************************
 * readDat is a prototype made to read in the LRLon, LRLat 
 * and zeroData files into memory
 * They are comma delimited files with doubles (lat and lon) and binary
 * All are representations of 1624 x 76 array of data
 * FDM 11-2013
 ************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX_ROWS 1624
#define MAX_WVC 76

char zeroName[200] = "/auto/users/minor/research/aveLatLon/c/zeroData.dat";
char latName[200] = "/auto/users/minor/research/aveLatLon/c/LRLat.dat";
char lonName[200] = "/auto/users/minor/research/aveLatLon/c/LRLon.dat";

void readInt(char* fileName, unsigned char data[][MAX_WVC]);
void readDbl(char* fileName, double data[][MAX_WVC]);

int main(int argc, char* argv[]) {
	double LRLat[MAX_ROWS][MAX_WVC];
	double LRLon[MAX_ROWS][MAX_WVC];
	unsigned char zeroData[MAX_ROWS][MAX_WVC];
	int row,cell;

	printf("Reading in zero data... ");	
	readInt(zeroName,zeroData);
	printf("done.\n");

	printf("Reading in lat data... ");
	readDbl(latName, LRLat);
	printf("done.\n");

	printf("Reading in lon data... ");
	readDbl(lonName, LRLon);
	printf("done.\n");


	/* Print out zeroData for testing*/
//	printf("Printing array: \n");
//	for(row = 0; row < MAX_ROWS - 1; row++) {
//		printf("%lf",LRLat[row][0]);
//		for(cell = 1; cell < MAX_WVC - 1; cell++) {
//			printf(",%f",LRLat[row][cell]);
//		}
//		printf("\n");
//	} 	

}

void readInt(char* fileName, unsigned char data[][MAX_WVC]){
	FILE* file;
	long fileSize;
	size_t res;
	int row,cell;
	char* buffer; /* buffer to hold file data */
	char *row_res, *cell_res; /* for tokenizing delimited data  */
	char *saveptr1, *saveptr2;

	/* Open filestream */	
	file = fopen(fileName,"r");
	if(file == NULL) {
		printf("Failed to open zero data file");
		exit(1);
	}

	/* Allocate space and read into buffer */
	fseek(file,0,SEEK_END);
	fileSize = ftell(file);	
	rewind(file);

	buffer = (char*) malloc(sizeof(char) * fileSize);
	if(buffer == NULL){
		printf("Memory Error \n");
		exit(1);
	}
	
	res = fread(buffer,1,fileSize,file);
	if(res != fileSize){
		printf("Reading Error");
		exit(1);
	}
	
	/* Read 0s and 1s into zeroData array */
	row = 0;
	row_res = strtok_r(buffer, "\n", &saveptr1);
	while(row_res != NULL) {
		cell = 0;
		cell_res = strtok_r(row_res, ",", &saveptr2);
		while(cell_res != NULL) {
			data[row][cell] = atoi(cell_res);
			cell_res = strtok_r(NULL, ",", &saveptr2);
			cell++;
		}
		row++;
		row_res = strtok_r(NULL, "\n", &saveptr1);
	}	

	fclose(file);
	free(buffer);

}

void readDbl(char* fileName, double data[][MAX_WVC]){
	FILE* file;
	long fileSize;
	size_t res;
	int row,cell;
	char* buffer; /* buffer to hold file data */
	char *row_res, *cell_res; /* for tokenizing delimited data  */
	char *saveptr1, *saveptr2;

	/* Open filestream */	
	file = fopen(fileName,"r");
	if(file == NULL) {
		printf("Failed to open zero data file");
		exit(1);
	}

	/* Allocate space and read into buffer */
	fseek(file,0,SEEK_END);
	fileSize = ftell(file);	
	rewind(file);

	buffer = (char*) malloc(sizeof(char) * fileSize);
	if(buffer == NULL){
		printf("Memory Error \n");
		exit(1);
	}
	
	res = fread(buffer,1,fileSize,file);
	if(res != fileSize){
		printf("Reading Error");
		exit(1);
	}
	
	/* Read doubles into data array */
	row = 0;
	row_res = strtok_r(buffer, "\n", &saveptr1);
	while(row_res != NULL) {
		cell = 0;
		cell_res = strtok_r(row_res, ",", &saveptr2);
		while(cell_res != NULL) {
			data[row][cell] = atof(cell_res);
			cell_res = strtok_r(NULL, ",", &saveptr2);
			cell++;
		}
		row++;
		row_res = strtok_r(NULL, "\n", &saveptr1);
	}	

	fclose(file);
	free(buffer);
}






