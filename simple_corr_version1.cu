#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
int64_t getnumlines(const char *fname,const char comment);
int64_t read_ascii_file(const char *filename,double **xpos,double **ypos,double **zpos);
int setup_bins_double(const char *fname,double *rmin,double *rmax,int *nbin,double **rupp);

int64_t getnumlines(const char *fname,const char comment)
{
  //FILE *fp = ((void *)0);
  FILE *fp = NULL;
  const int MAXLINESIZE = 10000;
  int64_t nlines = 0;
  char str_line[MAXLINESIZE];
  fp = fopen(fname,"rt");
  //if (fp == ((void *)0)) {
  if (fp == NULL) {
    fprintf(stderr,"Error: Could not open file `%s'\n",fname);
    //perror(((void *)0));
    perror(NULL);
    return (- 1);
  }
  while(1){
    //if (fgets(str_line,MAXLINESIZE,fp) != ((void *)0)) {
//WARNING: this does not remove white-space. You might
//want to implement that (was never an issue for me)
    if (fgets(str_line,MAXLINESIZE,fp) != (NULL)) {
      if (str_line[0] != comment) 
        nlines++;
    }
     else 
      break; 
  }
  fclose(fp);
  return nlines;
}

int64_t read_ascii_file(const char *filename,double **xpos,double **ypos,double **zpos)
{
  int64_t numlines = getnumlines(filename,'#');
  if (numlines <= 0) 
    return numlines;
    //Ritu: besides change (void *) to NULL, typecasting in the three lines below
  double *x = (double *) (calloc(numlines,sizeof(( *x))));
  double *y = (double *) (calloc(numlines,sizeof(( *y))));
  double *z = (double *) (calloc(numlines,sizeof(( *z))));
  //if (x == ((void *)0) || y == ((void *)0) || z == ((void *)0)) {
  if (x == NULL || y == NULL || z == NULL) {
    free(x);
    free(y);
    free(z);
    fprintf(stderr,"Error: Could not allocate memory for %ld elements for the (x/y/z) arrays\n",numlines);
    //perror(((void *)0));
    perror(NULL);
    return (- 1);
  }
  FILE *fp = fopen(filename,"rt");
  //if (fp == ((void *)0)) {
  if (fp == NULL) {
    fprintf(stderr,"Error:Could not open file `%s' in function %s\n",filename,__FUNCTION__);
    fprintf(stderr,"This is strange because the function `getnumlines' successfully counted the number of lines in that file\n");
    fprintf(stderr,"Did that file (`%s') just get deleted?\n",filename);
    //perror(((void *)0));
    perror(NULL);
    return (- 1);
  }
  int64_t index = 0;
  const int nitems = 3;
  const int MAXLINESIZE = 10000;
  char buf[MAXLINESIZE];
  while(1){
    if (fgets(buf,MAXLINESIZE,fp) != ((void *)0)) {
      int nread = sscanf(buf,"%lf %lf %lf",&x[index],&y[index],&z[index]);
      if (nread == nitems) {
        index++;
      }
    }
     else {
      break; 
    }
  }
  fclose(fp);
  if (index != numlines) {
    fprintf(stderr,"Error: There are supposed to be `%'ld lines of data in the file\n",numlines);
    fprintf(stderr,"But could only parse `%'ld lines containing (x y z) data\n",index);
    fprintf(stderr,"exiting...\n");
    return (- 1);
  }
   *xpos = x;
   *ypos = y;
   *zpos = z;
  return numlines;
}

int setup_bins_double(const char *fname,double *rmin,double *rmax,int *nbin,double **rupp)
{
//set up the bins according to the binned data file
//the form of the data file should be <rlow  rhigh ....>
  const int MAXBUFSIZE = 1000;
  char buf[MAXBUFSIZE];
  //FILE *fp = ((void *)0);
  FILE *fp = NULL;
  double low;
  double hi;
  const char comment = '#';
  const int nitems = 2;
  int nread = 0;
   *nbin = ((int )(getnumlines(fname,comment))) + 1;
  // *rupp = (calloc(( *nbin + 1),sizeof(double )));
  *rupp = (double *)(calloc(( *nbin + 1),sizeof(double )));
  if (rupp == ((void *)0)) {
    fprintf(stderr,"Error: Could not allocate memory for %d bins to store the histogram limits\n", *nbin + 1);
    //perror(((void *)0));
    perror(NULL);
    return 1;
  }
  fp = fopen(fname,"rt");
  if (fp == ((void *)0)) {
    free(( *rupp));
    fprintf(stderr,"Error: Could not open file `%s'..exiting\n",fname);
    //perror(((void *)0));
    perror(NULL);
    return 1;
  }
  int index = 1;
  while(1){
    if (fgets(buf,MAXBUFSIZE,fp) != ((void *)0)) {
      nread = sscanf(buf,"%lf %lf",&low,&hi);
      if (nread == nitems) {
        if (index == 1) {
           *rmin = low;
          ( *rupp)[0] = low;
        }
        ( *rupp)[index] = hi;
        index++;
      }
    }
     else {
      break; 
    }
  }
   *rmax = ( *rupp)[index - 1];
  fclose(fp);
  ( *rupp)[ *nbin] =  *rmax;
  ( *rupp)[ *nbin - 1] =  *rmax;
  return 0;
}

/*
__device__ double myAtomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
*/

void __global__ kernel0(double * xpos,double * ypos,double * zpos,double * rupp,int64_t Npart,const double sqr_rmin,const double sqr_rmax,const int logbins,const double log10rmin,const double inv_dlogr,int nbins, int device_M , int device_N, int64_t* totalNpairs,double* totalRpavg){

 int64_t i =  blockIdx.x * blockDim.x + threadIdx.x;
 //int64_t tid = threadIdx.x;

 if(i< Npart) {
  for(int64_t j = (i+1);j < Npart;j++) {
      const double dx = xpos[i] - xpos[j];
      const double dy = ypos[i] - ypos[j];
      const double dz = zpos[i] - zpos[j];
      const double r2 = dx * dx + dy * dy + dz * dz;
      if(r2 < sqr_rmin || r2 >= sqr_rmax) 
         continue;
      const double r = sqrt(r2);
      if(logbins) {
        const int kbin =(int )((log10(r) - log10rmin) * inv_dlogr);
        totalNpairs[(i*nbins)+kbin + 1]++;
        totalRpavg[(i*nbins)+kbin + 1] += r;
      }else {
        for(int kbin = nbins - 1;kbin >= 1;kbin--) {
            if(r >= rupp[kbin - 1]) {
                totalNpairs[(i*nbins)+kbin]++;
                totalRpavg[(i*nbins)+kbin] += r;
                break;
            }
        }
      }
  }

  //__syncthreads();

 }

}

int main(int argc,char **argv)
{
  double *device_rupp;
 // double *device_rpavg;
 // int64_t *device_npairs;
  double *device_zpos;
  double *device_ypos;
  double *device_xpos;
  int64_t *device_totalNpairs;
  double *device_totalRpavg;
  if (argc < 3) {
    fprintf(stderr,"\n\tUsage: %s `filename (string)' `filename-with-bins (string)' `[log bins (boolean)]'\n\n",argv[0]);
    fprintf(stderr,"\t************************************************************\n");
    fprintf(stderr,"\tRequired\n");
    fprintf(stderr,"\t--------\n");
    fprintf(stderr,"\t filename                string, an ascii file containing particle data (white-space-separated, 3 columns of x y z)\n");
    fprintf(stderr,"\t filename-with-bins      string, an ascii file containing <rlow rmax> specifying logarithmic bins (number of lines equal the number of bins)\n");
    fprintf(stderr,"\n\tOptional\n");
    fprintf(stderr,"\t--------\n");
    fprintf(stderr,"\t log-bins                boolean, default 0. Supply `1' indicating that the supplied bins are logarithmic (assumed to be log10)\n");
    fprintf(stderr,"\t************************************************************\n\n");
    return 1;
  }
  int logbins = 0;
  if (argc > 3) {
    logbins = atoi(argv[3]);
    fprintf(stderr,"Assuming that bins are logarithmic. Using logbins = %d\n",logbins);
  }
  double *xpos;
  double *ypos;
  double *zpos;
  int64_t Npart = read_ascii_file(argv[1],&xpos,&ypos,&zpos);
  if (Npart <= 0) {
    return Npart;
  }
  double rmin;
  double rmax;
  double *rupp;
  int nbins;
  int status = setup_bins_double(argv[2],&rmin,&rmax,&nbins,&rupp);
  if (status < 0) {
    return status;
  }
  //double *rpavg = (calloc(nbins,sizeof(( *rpavg))));
  //int64_t *npairs = (calloc(nbins,sizeof(( *npairs))));

  double *rpavg = (double *) (calloc(nbins,sizeof(( *rpavg))));
  int64_t *npairs = (int64_t *) (calloc(nbins,sizeof(( *npairs))));

  double *totalRpavg = (double *) (calloc((nbins*Npart),sizeof(( *rpavg))));
  int64_t *totalNpairs = (int64_t *) (calloc((nbins*Npart),sizeof(( *npairs))));

  const double sqr_rmin = rmin * rmin;
  const double sqr_rmax = rmax * rmax;
  const double log10rmin = log10(rmin);
  const double log10rmax = log10(rmax);
/* because of the way nbins is implemented
       bin `0' is underflow, and bin `nbin' is overflow  */
  const double dlogr = (log10rmax - log10rmin) / (nbins - 1);
  const double inv_dlogr = 1.0 / dlogr;

//Please note this is the section wherein the number of blocks and threads are calculated.  To change the number of threads alter the dimBlock whereas to change the number of blocks alter the dimGrid

int D_rows = (Npart > 1024 ) ? Npart/1024 : Npart;
int D_cols = (Npart > 1024 ) ? 1024 : 1;
//Ritu:updating the number of D_rows
if ( Npart % 1024){
 D_rows++;
}

//printf("\nD_rows:%d, D_cols:%d\n",D_rows, D_cols);

dim3 dimGrid(D_rows,1);
dim3 dimBlock(D_cols,1);
cudaMalloc((void **) &device_xpos,(Npart)*sizeof(int64_t));
cudaMemcpy(device_xpos,xpos,(Npart)*sizeof(double),cudaMemcpyHostToDevice);
cudaMalloc((void **) &device_ypos,(Npart)*sizeof(double));
cudaMemcpy(device_ypos,ypos,(Npart)*sizeof(double),cudaMemcpyHostToDevice);
cudaMalloc((void **) &device_zpos,(Npart)*sizeof(double));
cudaMemcpy(device_zpos,zpos,(Npart)*sizeof(double),cudaMemcpyHostToDevice);

//cudaMalloc((void **) &device_npairs,(nbins)*sizeof(int64_t));
//cudaMemcpy(device_npairs,npairs,(nbins)*sizeof(int64_t),cudaMemcpyHostToDevice);

//cudaMalloc((void **) &device_rpavg,(nbins)*sizeof(double));
//cudaMemcpy(device_rpavg,rpavg,(nbins)*sizeof(double),cudaMemcpyHostToDevice);

cudaMalloc((void **) &device_rupp,(nbins)*sizeof(double));
cudaMemcpy(device_rupp,rupp,(nbins)*sizeof(double),cudaMemcpyHostToDevice);

cudaMalloc((void **) &device_totalNpairs,(nbins*Npart)*sizeof(int64_t));

cudaMalloc((void **) &device_totalRpavg,(nbins*Npart)*sizeof(double));

for (int64_t j=0; j<nbins*Npart; j++){
    totalNpairs[j]  = 0;
    totalRpavg[j] = 0.0;
}

cudaMemcpy(device_totalNpairs,totalNpairs,(nbins*Npart)*sizeof(int64_t),cudaMemcpyHostToDevice);
cudaMemcpy(device_totalRpavg,totalNpairs,(nbins*Npart)*sizeof(double),cudaMemcpyHostToDevice);

kernel0<<<dimGrid,dimBlock>>>(device_xpos,device_ypos,device_zpos,device_rupp,Npart,sqr_rmin,sqr_rmax,logbins,log10rmin,inv_dlogr,nbins,1,Npart, device_totalNpairs,device_totalRpavg);

//kernel0<<<dimGrid,dimBlock>>>(device_xpos,device_ypos,device_zpos,device_npairs,device_rpavg,device_rupp,Npart,sqr_rmin,sqr_rmax,logbins,log10rmin,inv_dlogr,nbins,1,Npart, device_totalNpairs,device_totalRpavg);


//size_t sharedMemory = (Npart*nbins*sizeof(int64_t)) + (Npart*nbins*sizeof(double));

//kernel0<<<dimGrid,dimBlock, sharedMemory>>>(device_xpos,device_ypos,device_zpos,device_npairs,device_rpavg,device_rupp,Npart,sqr_rmin,sqr_rmax,logbins,log10rmin,inv_dlogr,nbins,totalNpairs,totalRpavg);
/*
  int IPT_function_replace;
*/

//cudaDeviceSynchronize();

cudaFree(device_xpos);
cudaFree(device_ypos);
cudaFree(device_zpos);
//cudaFree(device_npairs);
//cudaFree(device_rpavg);


 cudaMemcpy(totalNpairs,device_totalNpairs,(nbins*Npart)*sizeof(int64_t), cudaMemcpyDeviceToHost);
 cudaFree(device_totalNpairs);

 cudaMemcpy(totalRpavg,device_totalRpavg,(nbins*Npart)*sizeof(double), cudaMemcpyDeviceToHost);
 cudaFree(device_totalRpavg);

 for(int64_t k=0; k<Npart; k++){
   for (int64_t j=0; j<nbins; j++){
    npairs[j] += totalNpairs[(k*nbins)+ j];
    rpavg[j] += totalRpavg[(k*nbins)+ j];
   }
  }


  cudaFree(device_rupp);
  double rlow = rupp[0];
  for (int i = 1; i < nbins; i++) {
    if (npairs[i] > 0) 
      rpavg[i] /= npairs[i];
    fprintf(stdout,"%e\t%e\t%e\t%12lu\t%e\n",rlow,rupp[i],rpavg[i],npairs[i],0.0);
    rlow = rupp[i];
  }

  free(xpos);
  free(ypos);
  free(zpos);
  free(rupp);
  free(npairs);
  free(rpavg);
  free(totalNpairs);
  free(totalRpavg);
  return 0;
}
