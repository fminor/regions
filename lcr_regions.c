/*****************************************************************
 * Filename:  make_lcr_regions.c
 * read information from a qscat2b file 
 * (find land in L2B file to 
 * 	make regions to process using LCR)
 * Usage: make_lcr_regions <L2B filename> <Sector output filename>
 *
 * written by FDM at BYU 10/6/1013 based on read_qscat2b_2.c that was
 * written by DGL at BYU 1/20/2003 based on read_qscat2b.c
 *
 * code written to find land and make regions is marked by the 
 * initials FDM in the comments (at least until the end of main)
 * or F. D. Minor in function declaration comments
 *****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h> /* fmod */
#include <hdf.h>

#define MAX_ROWS 1624
#define MAX_WVC 76
#define MAX_AMBIG 4
#define T_LENGTH 23

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define rabs(x) ((x) >= 0 ? (x) : -(x))

int32 sd_id;

typedef struct{
  int rank,dim[3];
  float *data;
} qsdata;

typedef struct{
  int rank,dim[3];
  int *data;
} qsdatai;

/* function prototypes */
void get_timetags(char infile[1024],char timetags[MAX_ROWS][T_LENGTH]);
void *get_mem(int32 dt,int32 numval);
void getTextAfterNewline();

/* function prototypes for targeting near coastal regions FDM */
unsigned char coastline_distance(float, float);
void readInt(char* fileName, unsigned char data[][MAX_WVC]);
void readDbl(char* fileName, double data[][MAX_WVC]);
void transposeLon(float*, double[][MAX_WVC]);
void fill_array(double[][MAX_WVC],double[][MAX_WVC],
		unsigned char [][MAX_WVC], int [][MAX_WVC]);
void fill_array2(float*, float*, int*, int[][MAX_WVC]);
void print_array(int [][MAX_WVC]);

/* global variables for reading hdf file*/
qsdata extract_sds(char *in_var,int irec,int slab_size);
qsdatai extract_sds_i(char *in_var,int irec,int slab_size);

/* global variables for near-coastal threshold */
int Coast_distance_threshold=20;  /* default threshold in km */
char landname[200]="/auto/share/ref/WorldLandDistMap.dat";
int NSX2_coast=36000, NSY2_coast=18000;
int NS_coast=36000*18000; /* =NSX2_coast*NSY2_coast */
unsigned char *landdis;

/* global variables for lat/lon data files */
char zeroName[200] = "/auto/users/minor/src/data/zeroData.dat";
char latName[200] = "/auto/users/minor/src/data/LRLat.dat";
char lonName[200] = "/auto/users/minor/src/data/LRLon.dat";

int main(int argc,char *argv[])
{
  char *filename,fn[200];
  char *sectorname,sn[200];
  char line[200],attr_name[512],attr_data2[512],attr_chk[512];
  int32 retn,attr_index,adata_type,count;
  int8 *attr_data;
  int chk, chk1;

  char timetags[MAX_ROWS][T_LENGTH];
  FILE *imf;

  qsdatai wvc_row,wvc_index;
  qsdata wvc_lat,wvc_lon;
  qsdatai wvc_quality_flag;

  int ir,slab_size;
  int itmp,irec1,irec2,incr2,incr3,iamb;
  int i,j;
  int iwvc_quality_flag;

  /* Local variables for targeting near coastal regions FDM */
  double LRLat[MAX_ROWS][MAX_WVC]; /* Estimated latitude */ 
  double LRLon[MAX_ROWS][MAX_WVC]; /* Estimated longitude */
  unsigned char zeroData[MAX_ROWS][MAX_WVC]; /* 1's represent areas unfit for processing */
  int landmap[MAX_ROWS][MAX_WVC]; /* For L2B distances to land (values 0-255)*/

  /* Read the input filename */
  if (argv[1] == NULL){
    printf("Enter the name of the input file:");
    fgets(line,sizeof(line),stdin);
    sscanf(line,"%s",fn);
    filename=fn;
  }
  else{
    filename=argv[1];
  }	
  printf("\nFILENAME: %s \n",filename);

  /* Read output filename */
  if (argv[2] == NULL){
    printf("Enter the name of the sector file:");
    fgets(line,sizeof(line),stdin);
    sscanf(line,"%s",fn);
    sectorname=fn;
  }
  else{
    sectorname=argv[2];
  }	


  /* Read the timetag info contained in the HDF VDATA */
  get_timetags(filename,timetags);

  /* Open the HDF input file and initiate the SD interface */
  sd_id = SDstart(filename, DFACC_RDONLY);

  /* Make sure that the file is a QuikSCAT Level 2B file */
  printf("Verifying file is a QuikSCAT L2B file...\n");
  attr_index=SDfindattr(sd_id,"ShortName");
  retn=SDattrinfo(sd_id,attr_index,attr_name, &adata_type,&count);
  attr_data=(int8 *)malloc(count * DFKNTsize(adata_type));
  retn=SDreadattr(sd_id,attr_index,attr_data);
  getTextAfterNewline(attr_data,2,attr_data2);
  (void)strcpy(attr_chk,"QSCATL2B");
  chk=strcmp(attr_chk,attr_data2);
  (void)strcpy(attr_chk,"SWSL2B");
  chk1=strcmp(attr_chk,attr_data2);

  if (chk != 0 && chk1 != 0){
    printf("ERROR: Input file is not a QuikSCAT/SeaWinds Level 2B file\n");
    return(0);
  }
  
  /* Load landdis, the Land Distance file*/
  printf("Reading land coast distance file: (%s)\n", landname);
  landdis = (unsigned char *) malloc(sizeof(unsigned char)*NS_coast);    
  if (landdis == NULL) {
      fprintf(stdout,"*** ERROR allocating land coast distance array ***\n");
      exit(-1);
  }
  imf = fopen(landname,"r");
  if (imf == NULL) {
    fprintf(stdout,"*** ERROR: cannot open land distance file %s\n",landname);
    exit(-1);
  }
  if (fread(landdis, sizeof(unsigned char), NS_coast, imf) != NS_coast) {
     fprintf(stdout,"*** ERROR reading land coast distance file\n");
     exit(-1);
  }
  fclose(imf);

  /* Read in wvc_row, wvc_lat, wvc_lon, and wvc_index in their entirety */
  /* Read each SDS in its entirety.  For an example of reading the */
  /* QuikSCAT SDS data in slabs, please refer to read_qscat2a.c.   */
  printf("Reading in L2B data...\n");
  ir=0;

  wvc_row=extract_sds_i("wvc_row",ir,MAX_ROWS);
  wvc_lat=extract_sds("wvc_lat",ir,MAX_ROWS);
  wvc_lon=extract_sds("wvc_lon",ir,MAX_ROWS);
  wvc_index=extract_sds_i("wvc_index",ir,MAX_ROWS);
  wvc_quality_flag=extract_sds_i("wvc_quality_flag",ir,MAX_ROWS);

  /* Read in LRLat, LRLon and zeroData into memory FDM */
  printf("Reading in zero data(%s)...\n", zeroName);	
  readInt(zeroName,zeroData);
  printf("Reading in lat data(%s)...\n", latName);
  readDbl(latName, LRLat);
  printf("Reading in lon data(%s)...\n", lonName);
  readDbl(lonName, LRLon);
  transposeLon(wvc_lon.data, LRLon);
  printf("Transpose reference: %f\n", LRLon[407][36]);

  /* Make 1624 x 76 array of land distances for L2B file FDM */
  printf("Making Distance array \n");
  fill_array(LRLat, LRLon, zeroData, landmap);
  print_array(landmap);
  //printf("Flag array \n");
  //fill_array2(wvc_lat.data, wvc_lon.data, wvc_quality_flag.data, landmap);
  //print_array(landmap);

  /* Make regions for sectorname file FDM */

  /* Free up malloc'd memory */
  free(wvc_row.data);
  free(wvc_lat.data);
  free(wvc_lon.data);
  free(wvc_index.data);
  free(landdis);
  free(wvc_quality_flag.data);

  wvc_row.data=NULL;
  wvc_lat.data=NULL;
  wvc_lon.data=NULL;
  wvc_index.data=NULL;
  landdis=NULL;
  wvc_quality_flag.data=NULL;

  SDend(sd_id);
}

 /***************************************************************
  * print_array
  * prints a 1624 x 76 array of ints
  * param array is the array we are printing
  * FDM
  ***************************************************************/ 
void print_array(int landmap[][MAX_WVC]){
	int i,j;
	for(i=0; i < MAX_ROWS; i++) {
		j = 0;
		printf("%d",landmap[i][j]);
		for(j = 1; j < MAX_WVC; j++) {
			printf(",%d", landmap[i][j]);
		}
		printf("\n");
	}
}


/*****************************************************************
 * fill_array
 * Makes 1624 x 76 array of land distances for L2B file 
 * Also: excludes outside 60S - 60N range
 * param qsdata wvc_lat
 * param qsdata wvc_lon
 * param int landmap[][MAX_WVC] is the matrix to fill with data
 * 	10/2013 F. D. Minor
 ****************************************************************/
void fill_array(double LRLat[][MAX_WVC], double LRLon[][MAX_WVC], 
		unsigned char zeroData[][MAX_WVC], int landmap[][MAX_WVC]){
  float lat,lon;
  unsigned char dis;
  int ir,j;
  for(ir=0; ir < MAX_ROWS; ir++) {
	for(j=0; j < MAX_WVC; j++) {
		if(!zeroData[ir][j]) {
			lat = (float) LRLat[ir][j];
			lon = (float) LRLon[ir][j];
			dis = coastline_distance(lat,lon);
			landmap[ir][j] = dis;
		} else {
			landmap[ir][j] = -1;
		}		
	}
  }
}

/*****************************************************************
 * transposeLon
 * transposes our LR Longitude estimation to overlay the L2B swath
 * perfered reference is 407,36
 * range extends to 407,[24-46] (in first equitorial crossing
 * param qsdata wvc_lon
 * param double LRLon
 * 	11/2013 F. D. Minor
 ****************************************************************/
void transposeLon(float* wvc_lon, double LRLon[][MAX_WVC]){
	double lon,ref,temp;
	int ir,j,incr;
	ir = 407; /* First equitorial crossing */
	j = 36;
	incr = ir * MAX_WVC + j;
	ref = *(wvc_lon + incr);

	if(!ref) /* Find valid value as reference */
		for(j=24; j < 46; j++ && !ref) {
			incr = ir * MAX_WVC + j;
			ref = *(wvc_lon + incr);
		}
	if(!ref) {
		ir = 1218; /* Second equitorial crossing */ 
		j = 36;
		incr = ir * MAX_WVC + j;
		ref = *(wvc_lon + incr);
	}
	if(!ref) 
		for(j=24; j < 46; j++ && !ref) {
			incr = ir * MAX_WVC + j;
			ref = *(wvc_lon + incr);
		}
	printf("Reference longitude at (%d,%d) is: %f\n", ir,j,ref);
	
	/* Adjust reference longitude (if not from 407,36) */
	lon = LRLon[ir][j];
	printf("Here, normalized lon is: %f\n", lon);
	ref = ref - lon;
	printf("New reference longitude: %f\n", ref);

	/* Transpose LRLon by reference longitude */
	for(ir=0; ir < MAX_ROWS; ir++) {
		for(j=0; j < MAX_WVC; j++) {
			temp = LRLon[ir][j] + ref;
			if(temp > 360)
				temp = temp - 360;
			//printf("%d,", (int) temp);
			LRLon[ir][j] = temp;
		}
		//printf("\n");
	}

}


/*****************************************************************
 * fill_array2
 * Makes 1624 x 76 array of land flags for L2B file 
 * Also: excludes outside 60S - 60N range
 * param qsdata wvc_lat
 * param qsdata wvc_lon
 * param int landmap[][MAX_WVC] is the matrix to fill with data
 * 	10/2013 F. D. Minor
 ****************************************************************/
void fill_array2(float* wvc_lat, float* wvc_lon, int* wvc_quality, int landmap[][MAX_WVC]){
  float lat,lon;
  int land;
  int ir,j,incr;
  int mask = 1 << 7;  /* For isolating the 7th bit */
  printf("Entered fill_array2\n");
  for(ir=0; ir < MAX_ROWS; ir++) {
	for(j=0; j < MAX_WVC; j++) {
		//printf("R: %d C: %d ", ir, j);
		incr = ir*MAX_WVC+j;
		lat = *(wvc_lat+incr);
		lon = *(wvc_lon+incr);
		land = (*(wvc_quality+incr) & mask) ? 1 : 0;
		//land = *(wvc_quality+incr);
		//land = land & mask;
		//printf("Lat: %f Lon: %f ", lat, lon);
		//dist = coastline_distance(lat, lon);
		//printf("Flag: %x\n", land);
   		// !! Exclude values outside lat range
		landmap[ir][j] = land;		
	}
  }
}


/*****************************************************************
 * readInt
 * Reads in data from file and makes a 1624 x 76 array of ints
 * that in this case represent zeroData
 * the file is 1624 x 76 and is comma delimited 
 * param char* fileName - File containting data
 * param unsigned char data[][] - Array to store data
 * 	11/2013 F. D. Minor
 ****************************************************************/
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

/*****************************************************************
 * readDbl
 * Reads in data from file and makes a 1624 x 76 array of doubles
 * that in this case represent low resolution lat and lon
 * the file is 1624 x 76 and is comma delimited 
 * param char* fileName - File containting data
 * param double char data[][] - Array to store data
 * 	11/2013 F. D. Minor
 ****************************************************************/
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
		printf("\nFailed to open lat/lon data file\n");
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

/*****************************************************************
 * coastline_distance: 	a function that returns distance to land 
 * 			in km for any given lat and lon
 * param float: 	latitude
 * param float: 	longitude
 * return int: 		distance to land
 * 
 * 	Origionally written by D. Long
 * 	Copied to and adapted for this program:
 * 	10/2013	F.D. Minor
 ****************************************************************/	

unsigned char coastline_distance(float lat, float lon)
{
  int n,dis;
  int i,j;

  //if(((int) lat == 0) && ((int) lon == 0))
  //	return 0;

  i=(int) (fmod(lon+180.0,360.0)*100.0+0.5);
  j=(int) ((lat+90.0)*100.0+0.5);

  if (j<0) j=0;
  if (j>17999) j=17999;
  if (i<0) i=0;
  if (i>37999) i=0;

  n = j*NSX2_coast+i;
  dis = *(landdis+n);

  return(dis);
}



/*****************************************************************
  GET_TIMETAGS:  a function to read the timetags which are 
                 stored as VDATA.

  2/1999 K.L. Perry
*****************************************************************/

void get_timetags(char infile[1024],char timetags[MAX_ROWS][T_LENGTH])
{
  intn status_n;
  int32 status_32,file_id,vdata_id,vdata_ref;
  char vdata_name[1024],fieldname_list[1024];
  int i;

  /* Open the HDF file */
  file_id=Hopen(infile,DFACC_RDONLY,0);

  /* Initialize the VS interface. */
  status_n=Vstart(file_id);

  /* Get the reference number for the first vdata in the file. */
  vdata_ref=-1;
  vdata_ref=VSgetid(file_id,vdata_ref);
  
  /* Attach to the vdata for reading if it is found, otherwise 
     exit the program. */

  if (vdata_ref == 0) exit(0);
  vdata_id = VSattach (file_id, vdata_ref, "r"); 

  for (i=0;i<MAX_ROWS;i++){
    status_32=VSread(vdata_id,(uint8 *)timetags[i],1,FULL_INTERLACE);
  }

  /* Terminate access to the vdata and to the VS interface, then close 
     the HDF file. */
  status_32=VSdetach(vdata_id);
  status_n=Vend(file_id);
  status_32=Hclose(file_id);

  /* Display the read data as many records as the number of records 
     returned by VSread. */
  /*
    printf ("%s %s %s %s %s\n",timetags[0],timetags[1],timetags[1624/2-1],
             timetags[1624/2],timetags[1624-1]);
  */
}

/*****************************************************************
  EXTRACT_SDS:  a function to read the data which are stored
	        as Scientific Data Sets (SDS).

  2/1999 K.L. Perry, D.A. Seal
  1/2003 D.G. Long  + made more efficient for QuikSCAT/SeaWinds case

*****************************************************************/

qsdata extract_sds(char *in_var,int irec,int slab_size)
{
  char name[512];
  int32 sds_index,sds_id,retn;
  int32 rank,data_type,dim_sizes[3],num_type,nattrs;
  int32 start[3],edges[3];
  int32 i,numval;
  float64 cal,cal_err,off,off_err;
  void *out_var;
  qsdata caldata;
  int iflag=0;

  int8 *i8;
  uint8 *u8;
  int16 *i16;
  uint16 *u16;
  int32 *i32;
  uint32 *u32;
  float32 *f32;
  float64 *f64;
  float *ftmp;

  dim_sizes[0]=dim_sizes[1]=dim_sizes[2]=0;
  start[0]=irec;
  start[1]=start[2]=0;
  
  /* Search for the index of "in_var" */
  sds_index=SDnametoindex(sd_id,in_var);
  
  /* Select data set corresponding to the returned index */
  sds_id=SDselect(sd_id,sds_index);
  retn=SDgetinfo(sds_id,name,&rank,dim_sizes,&data_type,&nattrs);

  edges[0]=slab_size;
  edges[1]=dim_sizes[1];
  edges[2]=dim_sizes[2];
  
  caldata.rank=rank;
  caldata.dim[0]=dim_sizes[0];
  caldata.dim[1]=dim_sizes[1];
  caldata.dim[2]=dim_sizes[2];

  for (i=0,numval=1;i<rank;i++) numval *=caldata.dim[i];

  /* Allocate the memory for the output variable */
  out_var=get_mem(data_type,numval);

  /* Read the data set into the "out_var" array */
  retn=SDreaddata(sds_id,start,NULL,edges,(VOIDP)out_var);

  /* Get the calibration and offset values of input */
  retn=SDgetcal(sds_id,&cal,&cal_err,&off,&off_err,&num_type);

  /* Print Local Attributes */
  /*
    printf("%-19s %-3d %-4d %-4d %-4d %-2.9f %-2.9f %3d\n",
            name,rank,dim_sizes[0],dim_sizes[1],dim_sizes[2], cal,
            off,data_type);
   */

  retn=SDendaccess(sds_id);

  if (cal == 1.0) iflag=1;
  /*
  printf("in_var: %s  %d cal=%f %d  %d %d\n",in_var,data_type,cal,iflag,rank,numval);
  */

  /* Allocate memory for callibrated output */
  caldata.data=malloc(numval*sizeof(float));

  if (data_type == 20){
    if (iflag==1)
      for (i=0,i8=out_var,ftmp=caldata.data; i<numval;i++,i8++,ftmp++)
	*ftmp =*i8;
    else
      for (i=0,i8=out_var,ftmp=caldata.data; i<numval;i++,i8++,ftmp++)
	*ftmp =*i8*cal;
  }
  else if (data_type == 21){
    if (iflag==1)
      for (i=0,u8=out_var,ftmp=caldata.data; i<numval;i++,u8++,ftmp++)
	*ftmp =*u8;
    else
      for (i=0,u8=out_var,ftmp=caldata.data; i<numval;i++,u8++,ftmp++)
	*ftmp =*u8*cal;
  }
  else if (data_type == 22){
    if (iflag==1)
      for (i=0,i16=out_var,ftmp=caldata.data; i<numval;i++,i16++,ftmp++)
	*ftmp =*i16;
    else
      for (i=0,i16=out_var,ftmp=caldata.data; i<numval;i++,i16++,ftmp++)
	*ftmp =*i16*cal;
  }
  else if (data_type == 23){
    if (iflag==1)
      for (i=0,u16=out_var,ftmp=caldata.data; i<numval;i++,u16++,ftmp++)
	*ftmp =*u16;
    else
      for (i=0,u16=out_var,ftmp=caldata.data; i<numval;i++,u16++,ftmp++)
	*ftmp =*u16*cal;
  }
  else if (data_type == 24){
    if (iflag==1)
      for (i=0,i32=out_var,ftmp=caldata.data; i<numval;i++,i32++,ftmp++)
	*ftmp =*i32;
    else
      for (i=0,i32=out_var,ftmp=caldata.data; i<numval;i++,i32++,ftmp++)
	*ftmp =*i32*cal;
  }
  else if (data_type == 25){   /* potential conversion problem, but OK */
    if (iflag==1)
      for (i=0,u32=out_var,ftmp=caldata.data; i<numval;i++,u32++,ftmp++)
	*ftmp =*u32*cal;
    else
      for (i=0,u32=out_var,ftmp=caldata.data; i<numval;i++,u32++,ftmp++)
	*ftmp =*u32*cal;
  }
  else if (data_type == 5){
    if (iflag==1)
      for (i=0,f32=out_var,ftmp=caldata.data; i<numval;i++,f32++,ftmp++)
	*ftmp =*f32;
    else
      for (i=0,f32=out_var,ftmp=caldata.data; i<numval;i++,f32++,ftmp++)
	*ftmp =*f32*cal;
  }
  else if (data_type == 6){  /* poor conversion, but not used */
    if (iflag==1)
      for (i=0,f64=out_var,ftmp=caldata.data; i<numval;i++,f64++,ftmp++)
	*ftmp =*f64;
    else
      for (i=0,f64=out_var,ftmp=caldata.data; i<numval;i++,f64++,ftmp++)
	*ftmp =*f64*cal;
  }
  else{
    printf("Unknown Data Type \n");
  }

  return(caldata);
}

/*****************************************************************
  EXTRACT_SDS_I:  a function to read the data which are stored
	        as Scientific Data Sets (SDS).

  2/1999 K.L. Perry, D.A. Seal
  1/2003 D.G. Long  + made more efficient for QuikSCAT/SeaWinds case

*****************************************************************/

qsdatai extract_sds_i(char *in_var,int irec,int slab_size)
{
  char name[512];
  int32 sds_index,sds_id,retn;
  int32 rank,data_type,dim_sizes[3],num_type,nattrs;
  int32 start[3],edges[3];
  int32 i,numval;
  float64 cal,cal_err,off,off_err;
  void *out_var;
  qsdatai caldata;
  int iflag=0;

  int8 *i8;
  uint8 *u8;
  int16 *i16;
  uint16 *u16;
  int32 *i32;
  uint32 *u32;
  float32 *f32;
  float64 *f64;
  int *ftmp;

  dim_sizes[0]=dim_sizes[1]=dim_sizes[2]=0;
  start[0]=irec;
  start[1]=start[2]=0;
  
  /* Search for the index of "in_var" */
  sds_index=SDnametoindex(sd_id,in_var);
  
  /* Select data set corresponding to the returned index */
  sds_id=SDselect(sd_id,sds_index);
  retn=SDgetinfo(sds_id,name,&rank,dim_sizes,&data_type,&nattrs);

  edges[0]=slab_size;
  edges[1]=dim_sizes[1];
  edges[2]=dim_sizes[2];
  
  caldata.rank=rank;
  caldata.dim[0]=dim_sizes[0];
  caldata.dim[1]=dim_sizes[1];
  caldata.dim[2]=dim_sizes[2];

  for (i=0,numval=1;i<rank;i++) numval *=caldata.dim[i];

  /* Allocate the memory for the output variable */
  out_var=get_mem(data_type,numval);

  /* Read the data set into the "out_var" array */
  retn=SDreaddata(sds_id,start,NULL,edges,(VOIDP)out_var);

  /* Get the calibration and offset values of input */
  retn=SDgetcal(sds_id,&cal,&cal_err,&off,&off_err,&num_type);

  /* Print Local Attributes */
  /*
    printf("%-19s %-3d %-4d %-4d %-4d %-2.9f %-2.9f %3d\n",
            name,rank,dim_sizes[0],dim_sizes[1],dim_sizes[2], cal,
            off,data_type);
   */

  retn=SDendaccess(sds_id);

  if (cal == 1.0) iflag=1;
  /*
  printf("in_var: %s  %d cal=%f %d  %d %d\n",in_var,data_type,cal,iflag,rank,numval);
  */

  /* Allocate memory for callibrated output */
  caldata.data=malloc(numval*sizeof(float));

  if (data_type == 20){
    if (iflag==1)
      for (i=0,i8=out_var,ftmp=caldata.data; i<numval;i++,i8++,ftmp++)
	*ftmp =*i8;
    else
      for (i=0,i8=out_var,ftmp=caldata.data; i<numval;i++,i8++,ftmp++)
	*ftmp =*i8*cal;
  }
  else if (data_type == 21){
    if (iflag==1)
      for (i=0,u8=out_var,ftmp=caldata.data; i<numval;i++,u8++,ftmp++)
	*ftmp =*u8;
    else
      for (i=0,u8=out_var,ftmp=caldata.data; i<numval;i++,u8++,ftmp++)
	*ftmp =*u8*cal;
  }
  else if (data_type == 22){
    if (iflag==1)
      for (i=0,i16=out_var,ftmp=caldata.data; i<numval;i++,i16++,ftmp++)
	*ftmp =*i16;
    else
      for (i=0,i16=out_var,ftmp=caldata.data; i<numval;i++,i16++,ftmp++)
	*ftmp =*i16*cal;
  }
  else if (data_type == 23){
    if (iflag==1)
      for (i=0,u16=out_var,ftmp=caldata.data; i<numval;i++,u16++,ftmp++)
	*ftmp =*u16;
    else
      for (i=0,u16=out_var,ftmp=caldata.data; i<numval;i++,u16++,ftmp++)
	*ftmp =*u16*cal;
  }
  else if (data_type == 24){
    if (iflag==1)
      for (i=0,i32=out_var,ftmp=caldata.data; i<numval;i++,i32++,ftmp++)
	*ftmp =*i32;
    else
      for (i=0,i32=out_var,ftmp=caldata.data; i<numval;i++,i32++,ftmp++)
	*ftmp =*i32*cal;
  }
  else if (data_type == 25){   /* potential conversion problem, but OK */
    if (iflag==1)
      for (i=0,u32=out_var,ftmp=caldata.data; i<numval;i++,u32++,ftmp++)
	*ftmp =*u32*cal;
    else
      for (i=0,u32=out_var,ftmp=caldata.data; i<numval;i++,u32++,ftmp++)
	*ftmp =*u32*cal;
  }
  else if (data_type == 5){
    if (iflag==1)
      for (i=0,f32=out_var,ftmp=caldata.data; i<numval;i++,f32++,ftmp++)
	*ftmp =*f32;
    else
      for (i=0,f32=out_var,ftmp=caldata.data; i<numval;i++,f32++,ftmp++)
	*ftmp =*f32*cal;
  }
  else if (data_type == 6){  /* poor conversion, but not used */
    if (iflag==1)
      for (i=0,f64=out_var,ftmp=caldata.data; i<numval;i++,f64++,ftmp++)
	*ftmp =*f64;
    else
      for (i=0,f64=out_var,ftmp=caldata.data; i<numval;i++,f64++,ftmp++)
	*ftmp =*f64*cal;
  }
  else{
    printf("Unknown Data Type \n");
  }

  return(caldata);
}


/*****************************************************************
   GET_MEM:  a function which allocates the memory for an SDS
              using the data type.

   12/1998 K.L. Perry
           -adapted from get_type.c by A.V. Tran
*****************************************************************/

void *get_mem(int32 dt,int32 numval)
{
  switch(dt){

  case DFNT_INT8   : return(malloc(numval*sizeof(int8)));
  case DFNT_UINT8  : return(malloc(numval*sizeof(uint8)));
  case DFNT_INT16  : return(malloc(numval*sizeof(int16)));
  case DFNT_UINT16 : return(malloc(numval*sizeof(uint16)));
  case DFNT_INT32  : return(malloc(numval*sizeof(int32)));
  case DFNT_UINT32 : return(malloc(numval*sizeof(uint32)));
  case DFNT_FLOAT32: return(malloc(numval*sizeof(float32)));
  case DFNT_FLOAT64: return(malloc(numval*sizeof(float64)));

  default : return("unknown number type");
  }
}

/*****************************************************************
   getTextAfterNewline: a function which extracts all characters 
                        between two new lines from a string.

   7/1998 D.A. Seal
*****************************************************************/

void getTextAfterNewline( char *inputtext, int nn, char outputtext[240] )
{
  int nni, ii;
  char *ch;
  
  for ( nni = 0, ii = 0, ch = inputtext; *ch != '\0'; ch++ ) {
    if ( *ch == '\n' ) nni++;
    if ( nni == nn && *ch != '\n' ) {
      outputtext[ii] = *ch;
      ii++;
    }
    if ( nni > nn ) break;
  }
  outputtext[ii] = '\0';
  return;
}
