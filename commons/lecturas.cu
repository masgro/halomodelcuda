__constant__ float  d_coef_nu_m[7];
__constant__ int  d_coef_nu_m_n;
__constant__ float  d_coef_dlognu_dlogm[6];
__constant__ int  d_coef_dlognu_dlogm_n;
__constant__ float  d_coef_m_nu[7];
__constant__ int  d_coef_m_nu_n;
__constant__ float  d_coef_xilin[7];
__constant__ int  d_coef_xilin_n;

float *h_coef_nu_m;
int    h_coef_nu_m_n;
float *h_coef_dlognu_dlogm;
int    h_coef_dlognu_dlogm_n;
float *h_coef_m_nu;
int    h_coef_m_nu_n;
float *h_coef_xilin;
int    h_coef_xilin_n;


__inline__ void lee_ajustes_alineamiento(void)
{
  int ii, nfile = 10;
  int j, nlines = 30;
  FILE *pfin;
  char filename[200];
  float norm;

  /*Define la masa minima y maxima para el alineamiento*/
  h_alig_m[0] = 12.82;
  h_alig_m[1] = 12.87;
  h_alig_m[2] = 12.92;
  h_alig_m[3] = 12.98;
  h_alig_m[4] = 13.05;
  h_alig_m[5] = 13.13;
  h_alig_m[6] = 13.22;
  h_alig_m[7] = 13.33;
  h_alig_m[8] = 13.49;
  h_alig_m[9] = 13.74;
  h_alig_m[10] = 15.11;

  /*Define el radio minimo y maximo para el alineamiento*/
  h_alig_rmin = 0.0f;
  h_alig_rmax = 2.0f;
  h_alig_dr   = (h_alig_rmax - h_alig_rmin)/(float)nlines;

  /*Lee los coeficientes del alineamiento*/
  for(ii=0;ii<nfile;ii++)
  {
    sprintf(filename,"align/ajuste_%02d.dat",ii);
    if(!(pfin=fopen(filename,"r")))
      {printf("no se pudo abrir el archivo %s\n",filename);exit(0);}
    for(j=0;j<nlines;j++)
    {
      if(!fscanf(pfin,"%f %f %f %f %f\n",
          &h_alig[ii][j][0],&h_alig[ii][j][1],&h_alig[ii][j][2],
          &h_alig[ii][j][3],&h_alig[ii][j][4]))
        {printf("no se pudo leer el archivo %s\n",filename);exit(0);}

      norm = h_alig[ii][j][0] + h_alig[ii][j][1]/2.0 + h_alig[ii][j][2]/3.0 + h_alig[ii][j][3]/4.0
             + h_alig[ii][j][4]/5.0;

      h_alig[ii][j][0] /= norm;
      h_alig[ii][j][1] /= norm; 
      h_alig[ii][j][2] /= norm;
      h_alig[ii][j][3] /= norm;
      h_alig[ii][j][4] /= norm;

    }
  }

  /*Copia los ajustes del alineamiento al device*/
  HANDLE_ERROR(cudaMemcpyToSymbol(d_alig,h_alig,10*30*5*sizeof(float)));

  HANDLE_ERROR(cudaMemcpyToSymbol(d_alig_m,h_alig_m,11*sizeof(float)));

  HANDLE_ERROR(cudaMemcpyToSymbol(d_alig_rmin,&h_alig_rmin,sizeof(float)));
  HANDLE_ERROR(cudaMemcpyToSymbol(d_alig_rmax,&h_alig_rmax,sizeof(float)));
  HANDLE_ERROR(cudaMemcpyToSymbol(d_alig_dr,&h_alig_dr,sizeof(float)));
}

/*Lee los parametros de los ajustes para la forma*/
__inline__ void lee_parametros(void)
{
  FILE *pfin;
  char filename[200];

  /*Lee los ajustes de la forma: parametro bc*/
  sprintf(filename,"formas/ajustebc.dat");
  if(!(pfin=fopen(filename,"r")))
    {printf("no se pudo abrir el archivo %s\n",filename);exit(0);}
  if(!fscanf(pfin,"%f %f %f\n",&h_bc[0][0],&h_bc[0][1],&h_bc[0][2]))
    {printf("no se pudo leer el archivo %s\n",filename);exit(0);}
  if(!fscanf(pfin,"%f %f %f\n",&h_bc[1][0],&h_bc[1][1],&h_bc[1][2]))
    {printf("no se pudo leer el archivo %s\n",filename);exit(0);}
  if(!fscanf(pfin,"%f %f %f\n",&h_bc[2][0],&h_bc[2][1],&h_bc[2][2]))
    {printf("no se pudo leer el archivo %s\n",filename);exit(0);}
  if(!fscanf(pfin,"%f %f %f\n",&h_bc[3][0],&h_bc[3][1],&h_bc[3][2]))
    {printf("no se pudo leer el archivo %s\n",filename);exit(0);}
  fclose(pfin);

  /*Lee los ajustes de la forma: parametro ab*/
  sprintf(filename,"formas/ajusteab.dat");
  if(!(pfin=fopen(filename,"r")))
    {printf("no se pudo abrir el archivo %s\n",filename);exit(0);}
  if(!fscanf(pfin,"%f %f %f\n",&h_ab[0][0],&h_ab[0][1],&h_ab[0][2]))
    {printf("no se pudo leer el archivo %s\n",filename);exit(0);}
  if(!fscanf(pfin,"%f %f %f\n",&h_ab[1][0],&h_ab[1][1],&h_ab[1][2]))
    {printf("no se pudo leer el archivo %s\n",filename);exit(0);}
  if(!fscanf(pfin,"%f %f %f\n",&h_ab[2][0],&h_ab[2][1],&h_ab[2][2]))
    {printf("no se pudo leer el archivo %s\n",filename);exit(0);}
  if(!fscanf(pfin,"%f %f %f\n",&h_ab[3][0],&h_ab[3][1],&h_ab[3][2]))
    {printf("no se pudo leer el archivo %s\n",filename);exit(0);}
  fclose(pfin);

  /*Lee los ajustes de la forma: normalizacion*/
  sprintf(filename,"formas/ajuste_norma.dat");
  if(!(pfin=fopen(filename,"r")))
    {printf("no se pudo abrir el archivo %s\n",filename);exit(0);}
  if(!fscanf(pfin,"%f %f %f %f\n",&h_norm[0],&h_norm[1],&h_norm[2],&h_norm[3]))
    {printf("no se pudo leer el archivo %s\n",filename);exit(0);}
  fclose(pfin);

  /*Copia los ajustes de la forma al device*/
  HANDLE_ERROR(cudaMemcpyToSymbol(d_bc,h_bc,4*3*sizeof(float)));
  HANDLE_ERROR(cudaMemcpyToSymbol(d_ab,h_ab,4*3*sizeof(float)));
  HANDLE_ERROR(cudaMemcpyToSymbol(d_norm,h_norm,4*sizeof(float)));
}

void read_coefficients(void)
{
  int  i, n;
  char filename[200];
  FILE *pfin;

  sprintf(filename,"coeficientes.dat");
  pfin = fopen(filename,"r");
  assert(pfin != NULL);

  fprintf(stdout,"Leyendo coeficientes...\n");

  /*Lee coeficientes Nu(Mass)*/
  fscanf(pfin,"%d\n",&n);
  h_coef_nu_m_n = n;
  h_coef_nu_m = (float *) malloc(n*sizeof(float));
  for(i = 0; i < n; i++)
    fscanf(pfin,"%f\n",&h_coef_nu_m[i]);
  HANDLE_ERROR(cudaMemcpyToSymbol(d_coef_nu_m_n,&h_coef_nu_m_n,sizeof(int)));
  HANDLE_ERROR(cudaMemcpyToSymbol(d_coef_nu_m,h_coef_nu_m,h_coef_nu_m_n*sizeof(float)));
  HANDLE_ERROR(cudaMemcpyFromSymbol(&i,d_coef_nu_m_n,sizeof(int)));
  assert(i == h_coef_nu_m_n);

  /*Fija coeficientes dlogNu_dlogMass*/
  h_coef_dlognu_dlogm_n = h_coef_nu_m_n - 1;
  h_coef_dlognu_dlogm = (float *) malloc(h_coef_dlognu_dlogm_n*sizeof(float));
  for(i = 0; i < h_coef_dlognu_dlogm_n; i++)
    h_coef_dlognu_dlogm[i] = (float)(i+1)*h_coef_nu_m[i+1];
  HANDLE_ERROR(cudaMemcpyToSymbol(d_coef_dlognu_dlogm_n,&h_coef_dlognu_dlogm_n,sizeof(int)));
  HANDLE_ERROR(cudaMemcpyToSymbol(d_coef_dlognu_dlogm,h_coef_dlognu_dlogm,h_coef_dlognu_dlogm_n*sizeof(float)));
  HANDLE_ERROR(cudaMemcpyFromSymbol(&i,d_coef_dlognu_dlogm_n,sizeof(int)));
  assert(i == h_coef_dlognu_dlogm_n);

  /*Lee coeficientes Mass(Nu)*/
  fscanf(pfin,"%d\n",&n);
  h_coef_m_nu_n = n;
  h_coef_m_nu = (float *) malloc(n*sizeof(float));
  for(i = 0; i < n; i++)
    fscanf(pfin,"%f\n",&h_coef_m_nu[i]);
  HANDLE_ERROR(cudaMemcpyToSymbol(d_coef_m_nu_n,&h_coef_m_nu_n,sizeof(h_coef_m_nu_n)));
  HANDLE_ERROR(cudaMemcpyToSymbol(d_coef_m_nu,h_coef_m_nu,h_coef_m_nu_n*sizeof(float)));
  HANDLE_ERROR(cudaMemcpyFromSymbol(&i,d_coef_m_nu_n,sizeof(int)));
  assert(i == h_coef_m_nu_n);

  /*Lee coeficientes Xilin*/
  fscanf(pfin,"%d\n",&n);
  h_coef_xilin_n = n;
  h_coef_xilin = (float *) malloc(n*sizeof(float));
  for(i = 0; i < n; i++)
    fscanf(pfin,"%f\n",&h_coef_xilin[i]);
  HANDLE_ERROR(cudaMemcpyToSymbol(d_coef_xilin_n,&h_coef_xilin_n,sizeof(h_coef_xilin_n)));
  HANDLE_ERROR(cudaMemcpyToSymbol(d_coef_xilin,h_coef_xilin,h_coef_xilin_n*sizeof(float)));
  HANDLE_ERROR(cudaMemcpyFromSymbol(&i,d_coef_xilin_n,sizeof(int)));
  assert(i == h_coef_xilin_n);
}

