/*****************************************************************
  Filename:  make_lcr_regions.c

  read information from a qscat2b file 
  (find land in L2B file to 
	make regions to process using LCR)
  Usage: make_lcr_regions <L2B filename> <Sector output filename>
  
  written by FDM at BYU 10/6/1013 based on read_qscat2b_2.c that was
  written by DGL at BYU 1/20/2003 based on read_qscat2b.c

******************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <hdf.h>

#define MAX_ROWS 1624
#define MAX_WVC 76
#define MAX_AMBIG 4
#define T_LENGTH 23

#define MAX_NORTH 60 
#define MAX_SOUTH 60

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
int coastline_distance(float, float);
void fill_array(float*, float*, int[][]);
void fill_array2(float*, float*, int*, int[][]);
void print_array(int[][]);

/* global variables for reading hdf file*/
qsdata extract_sds(char *in_var,int irec,int slab_size);
qsdatai extract_sds_i(char *in_var,int irec,int slab_size);

/* global variables for near-coastal threshold */
int Coast_distance_threshold=20;  /* default threshold in km */
char landname[200]="/auto/share/ref/WorldLandDistMap.dat";
int NSX2_coast=36000, NSY2_coast=18000;
int NS_coast=36000*18000; /* =NSX2_coast*NSY2_coast */
unsigned char *landdis;


int main(int argc,char *argv[])
{
  char *filename,fn[200];
  char *sectorname,sn[200];
  char line[200],attr_name[512],attr_data2[512],attr_chk[512];
  int32 retn,attr_index,adata_type,count;
  int8 *attr_data;
  int chk, chk1;

  char timetags[MAX_ROWS][T_LENGTH];
  int landmap[MAX_ROWS][MAX_WVC]; /* For L2B distances to land */
  FILE *imf;

  qsdatai wvc_row,wvc_index;
  qsdata wvc_lat,wvc_lon;
//  qsdatai num_in_fore,num_in_aft,num_out_fore,num_out_aft;
  qsdatai wvc_quality_flag;
//  qsdatai num_ambigs;
//  qsdata model_speed,model_dir,atten_corr;
//  qsdata wind_speed,wind_dir,wind_speed_err,wind_dir_err;
//  qsdata max_likelihood_est;
//  qsdatai wvc_selection;
//  qsdata wind_speed_selection,wind_dir_selection;
//  qsdata mp_rain_probability,nof_rain_index;

  int ir,slab_size;
  int itmp,irec1,irec2,incr2,incr3,iamb;
  int i,j;
  int iwvc_quality_flag;

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
    printf("Enter the name of the output file:");
    fgets(line,sizeof(line),stdin);
    sscanf(line,"%s",fn);
    sectorname=fn;
  }
  else{
    sectorname=argv[2];
  }	
  printf("\nSector File: %s \n",sectorname);


  /* Read the timetag info contained in the HDF VDATA */
  get_timetags(filename,timetags);

  /* Open the HDF input file and initiate the SD interface */
  sd_id = SDstart(filename, DFACC_RDONLY);

  /* Make sure that the file is a QuikSCAT Level 2B file */
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
  landdis = (unsigned char *) malloc(sizeof(unsigned char)*NS_coast);    
  if (landdis == NULL) {
      fprintf(stdout,"*** ERROR allocating land coast distance array ***\n");
      exit(-1);
  }
  printf("Open land coast distance file: %s\n",landname);
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
  ir=0;

  wvc_row=extract_sds_i("wvc_row",ir,MAX_ROWS);
  wvc_lat=extract_sds("wvc_lat",ir,MAX_ROWS);
  wvc_lon=extract_sds("wvc_lon",ir,MAX_ROWS);
  wvc_index=extract_sds_i("wvc_index",ir,MAX_ROWS);
//  num_in_fore=extract_sds_i("num_in_fore",ir,MAX_ROWS);
//  num_in_aft=extract_sds_i("num_in_aft",ir,MAX_ROWS);
//  num_out_fore=extract_sds_i("num_out_fore",ir,MAX_ROWS);
//  num_out_aft=extract_sds_i("num_out_aft",ir,MAX_ROWS);
  wvc_quality_flag=extract_sds_i("wvc_quality_flag",ir,MAX_ROWS);
//  atten_corr=extract_sds("atten_corr",ir,MAX_ROWS);
//  model_speed=extract_sds("model_speed",ir,MAX_ROWS);
//  model_dir=extract_sds("model_dir",ir,MAX_ROWS);
//  num_ambigs=extract_sds_i("num_ambigs",ir,MAX_ROWS);
//  wind_speed=extract_sds("wind_speed",ir,MAX_ROWS);
//  wind_dir=extract_sds("wind_dir",ir,MAX_ROWS);
//  wind_speed_err=extract_sds("wind_speed_err",ir,MAX_ROWS);
//  wind_dir_err=extract_sds("wind_dir_err",ir,MAX_ROWS);
//  max_likelihood_est=extract_sds("max_likelihood_est",ir,MAX_ROWS);
//  wvc_selection=extract_sds_i("wvc_selection",ir,MAX_ROWS);
//  wind_speed_selection=extract_sds("wind_speed_selection",ir,MAX_ROWS);
//  wind_dir_selection=extract_sds("wind_dir_selection",ir,MAX_ROWS);
//  mp_rain_probability=extract_sds("mp_rain_probability",ir,MAX_ROWS);
//  nof_rain_index=extract_sds("nof_rain_index",ir,MAX_ROWS);
  
  /* Select wind vector cell rows to be written to the screen */
/*  printf("\nEnter the first record number:");
  fgets(line,sizeof(line),stdin);
  sscanf(line,"%d",&irec1);

  printf("Enter the last record number:");
  fgets(line,sizeof(line),stdin);
  sscanf(line,"%d",&irec2);

  if (irec1 > irec2){
    itmp=irec1;
    irec1=irec2;
    irec2=itmp;
  }
  if ((irec1 < 1) || (irec2 > MAX_ROWS)){
    printf("ERROR: wvc rows must be between 1 and 1624 \n");
    return(0);
  }*/

  /* Subtract 1 from irec1 and irec2 to adjust for C running
     from 0 instead of 1 (so wvc_row 1 matches input of 1)   */

//  for(ir=irec1-1; ir<irec2; ir++){

    /* Print results to screen */

//    printf("\nTIME: %s\n",timetags[ir]);
//    printf("WVC ROW:  %5.0f \n",*(wvc_row.data+ir));

//    printf("WVC#  WVC_Qual WVC_Latitude/Longitude  Selected Wind Vector NWP Wind Vector Num/Sel Ambig  DRE Wind Vector  MP_Rain NOF_Rain\n");

/*    for (j=0; j<MAX_WVC; j++){
      incr2=ir*MAX_WVC+j;
      iamb=(*(wvc_selection.data+incr2));
      incr3=(ir*MAX_WVC*MAX_AMBIG)+(j*MAX_AMBIG)+(iamb-1);

      if(*(num_ambigs.data+incr2) > 0){
	iwvc_quality_flag=*(wvc_quality_flag.data+incr2);
	
	printf("%2.0d     0X%4.4x       %6.2f   %6.2f      %6.2f   %6.2f    %6.2f   %6.2f     %2.0d   %2.0d   %6.2f   %6.2f  %6.2f  %6.2f\n",
	       *(wvc_index.data+incr2),iwvc_quality_flag,
	       *(wvc_lat.data+incr2),*(wvc_lon.data+incr2),
	       *(wind_speed.data+incr3),*(wind_dir.data+incr3),
	       *(model_speed.data+incr2),*(model_dir.data+incr2),
	       *(num_ambigs.data+incr2),*(wvc_selection.data+incr2),
	       *(wind_speed_selection.data+incr2),
	       *(wind_dir_selection.data+incr2),
	       *(mp_rain_probability.data+incr2),
	       *(nof_rain_index.data+incr2));
      }
    }	*/
//  }

  /* Make 1624 x 76 array of land distances for L2B file */
  //printf("Distance array \n");
  //fill_array(wvc_lat.data, wvc_lon.data, landmap);
  //print_array(landmap);
  printf("Flag array \n");
  fill_array2(wvc_lat.data, wvc_lon.data, wvc_quality_flag.data, landmap);
  print_array(landmap);
  /* 
   * Make regions for sectorname file
   * 
   * */



  /* Free up malloc'd memory */
  free(wvc_row.data);
  free(wvc_lat.data);
  free(wvc_lon.data);
  free(wvc_index.data);
  free(landdis);
//  free(num_in_fore.data);
//  free(num_in_aft.data);
//  free(num_out_fore.data);
//  free(num_out_aft.data);
  free(wvc_quality_flag.data);
//  free(atten_corr.data);
//  free(model_speed.data);
//  free(model_dir.data);
//  free(num_ambigs.data);
//  free(wind_speed.data);
//  free(wind_dir.data);
//  free(wind_speed_err.data);
//  free(wind_dir_err.data);
//  free(max_likelihood_est.data);
//  free(wvc_selection.data);
//  free(wind_speed_selection.data);
//  free(wind_dir_selection.data);
//  free(mp_rain_probability.data);
//  free(nof_rain_index.data);

  wvc_row.data=NULL;
  wvc_lat.data=NULL;
  wvc_lon.data=NULL;
  wvc_index.data=NULL;
  landdis=NULL;
//  num_in_fore.data=NULL;
//  num_in_aft.data=NULL;
//  num_out_fore.data=NULL;
//  num_out_aft.data=NULL;
  wvc_quality_flag.data=NULL;
//  atten_corr.data=NULL;
//  model_speed.data=NULL;
//  model_dir.data=NULL;
//  num_ambigs.data=NULL;
//  wind_speed.data=NULL;
//  wind_dir.data=NULL;
//  wind_speed_err.data=NULL;
//  wind_dir_err.data=NULL;
//  max_likelihood_est.data=NULL;
//  wvc_selection.data=NULL;
//  wind_speed_selection.data=NULL;
//  wind_dir_selection.data=NULL;
//  mp_rain_probability.data=NULL;
//  nof_rain_index.data=NULL;

  SDend(sd_id);
}

 /***************************************************************
  * print_array
  * prints a 1624 x 76 array of ints
  * param array is the array we are printing
  *
  ***************************************************************/ 
void print_array(int landmap[][MAX_WVC]){
	int i,j;
	printf("Entered print array\n");
	for(i=0; i < MAX_ROWS; i++) {
		//printf("R: %d ", i);
		for(j = 0; j < MAX_WVC; j++) {
			printf("%d ", landmap[i][j]);
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
void fill_array(float* wvc_lat, float* wvc_lon, int landmap[][MAX_WVC]){
  float lat,lon;
  int dist;
  int ir,j,incr;
  printf("Entered fill_array\n");
  for(ir=0; ir < MAX_ROWS; ir++) {
	for(j=0; j < MAX_WVC; j++) {
		//printf("R: %d C: %d ", ir, j);
		incr = ir*MAX_WVC+j;
		lat = *(wvc_lat+incr);
		lon = *(wvc_lon+incr);
		//printf("Lat: %f Lon: %f ", lat, lon);
		dist = coastline_distance(lat, lon);
		//printf("Dist: %d\n", dist);
   		// !! Exclude values outside lat range
		landmap[ir][j] = dist;		
	}
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
  printf("Entered fill_array\n");
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

int coastline_distance(float lat, float lon)
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
  dis = (int) *(landdis+n);

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
