
/*

denseCOLOR compute histograms of color projection on a regular dense grid

Usage
------

[dcolor , infodcolor] = denseCOLOR(I , [options] );


Inputs
-------

I                                     Input image (ny x nx x [3]) in UINT8 format. 

options
	   scale                          Scaling vector (1 x nscale). Extract descriptors at different scaling of the image (default scale = [1]).
	   sigma_scale                    Scaling factor to obtain the standard deviation of the Gaussian filter (sigma = sigma_scale/scale)(default sigma_scale = 0.6)
	   deltax                         Division step in the x-axis for the grid (default deltax = floor(nx*min(scale))) 
	   deltay                         Division step in the y-axis for the grid (default deltay = floor(ny*min(scale)))
       color                          0 : force gray-scale (dimcolor = 1, default), 1 : RGB (dimcolor = 3), 2 : nRGB (dimcolor = 3), 3 : Opponent (dimcolor = 3), 
                                      4 : nOpponent (dimcolor = 2), 5 : Hue (dimcolor = 1)
	   nbins                          Number of bins for histograms (default nbins = 32)
	   patchsize                      Size of the patch where the descriptor is computed (default patchsize = nbins/2 )	  
	   norm                           Normalization : norm = 0 <=> no normalization, norm = 1 <=> v=v/(sum(v)+epsi), norm = 2 <=> v=v/sqrt(sum(v²)+epsi²), 
	                                  norm = 3 <=> v=sqrt(v/(sum(v)+epsi)) , norm = 3 <=> L2-clamped (default norm = 1)
	   clamp                          Clamping value (default clamp = 0.2)


Outputs
-------

dcolor                                COLOR descriptors (nbins x nb_pts) where nb_pts = deltax*deltay*nscale*dimcolor
infodcolor                            COLOR descriptors informations(7 x nb_pts)   where nb_pts = deltax*deltay*nscale*dimcolor
                                      infodcolor(1,i) = y
									  infodcolor(2,i) = x
									  infodcolor(3,i) = scale
									  infodcolor(4,i) = color
									  infodcolor(5,i) = nyscale
									  infodcolor(6,i) = nxscale
									  infodcolor(7,i) = currentpatchsize
									  


Example 1
---------


I                                    = imread('02769_Right_StudentOffice.jpeg');


options.scale                        = [0.5 , 0.75 , 1];
options.sigma_scale                  = 0.6;
options.deltax                       = 10;
options.deltay                       = 10;
options.color                        = 3;
options.nbins                        = 64;
options.patchsize                    = 24;
options.norm                         = 4;
options.clamp                        = 0.2;


[dcolor , infodcolor]                 = denseCOLOR(I , options ); 

figure(1)
imagesc(dcolor)


figure(2)
imagesc(I)
colormap(gray)
hold on
plot(infodcolor(2 , :) , infodcolor(1 , :) , 'r+')
hold off




Example 2
---------


I                                    = imread('image_0174.jpg');


options.scale                        = [1];
options.sigma_scale                  = 0.6;
options.deltax                       = 10;
options.deltay                       = 10;
options.nbins                        = 16;
options.patchsize                    = 16;
options.norm                         = 1;


[dcolor , infodcolor]                = denseCOLOR(I , options ); 

half                                 = options.patchsize/2;

figure(1)
imagesc(dcolor)

xr                                   = [infodcolor(2, :)-half ; infodcolor(2, :)-half ; infodcolor(2, :)+ half ; infodcolor(2, :)+ half ; infodcolor(2, :)-half] + 1.5;
yr                                   = [infodcolor(1, :)-half ; infodcolor(1, :)+half ; infodcolor(1, :)+ half ; infodcolor(1, :)- half ; infodcolor(1, :)-half] + 1.5;


figure(2)
imagesc(I)
colormap(gray)
hold on
plot(infodcolor(2 , :)+1.5 , infodcolor(1 , :)+1.5 , 'r+')
plot(xr , yr , 'b')
hold off


To compile
----------

mex  -output denseCOLOR.dll denseCOLOR.c

mex  -g -output denseCOLOR.dll denseCOLOR.c

mex  -f mexopts_intel10.bat -output denseCOLOR.dll denseCOLOR.c


Author : Sébastien PARIS : sebastien.paris@lsis.org
-------  Date : 02/04/2010


References :  
---------        


*/


#include <time.h>
#include <math.h>
#include <mex.h>


#define tiny 1e-8
#define verytiny 1e-15
#define PI 3.14159265358979323846
#define PIhalf2 1.5707963267948966
#define sqrt3 1.73205080756888
#define invsqrt2 0.707106781186547
#define invsqrt3 0.577350269189626
#define invsqrt6 0.408248290463863

#ifndef max
    #define max(a,b) (a >= b ? a : b)
    #define min(a,b) (a <= b ? a : b)
#endif
#define CTE 3.71692218884984   /* sqrt( log(10.0) * 6.0 ) )*/


struct opts
{

	double  *scale;
	int     nscale;
    double  sigma_scale;

	int    patchsize;
	int    deltax;
	int    deltay;

	int    color;
	int    nbins;

	int    norm;
	double clamp;
};


/*----------------------------------------------------------------------------------------------------------------------------------------- */
/*----------------------------------------------------------------------------------------------------------------------------------------- */
/*----------------------------------------------------------------------------------------------------------------------------------------- */

/* Function prototypes */

void rgb2gray(unsigned char * , int , int , double *);
void rgb2nrgb(unsigned char * , int , int , double *);
void rgb2opponent(unsigned char * , int , int , double *);
void rgb2nopponent(unsigned char * , int , int , double *);
void rgb2hue(unsigned char * , int , int , double *);

void gaussian_kernel( double *, double  , int  , int  );
void gaussian_sampler( double *  , int , int , double , double , double * , double *, double *);
void denseCOLOR(double * , int  , int , struct opts , double * , double *  , int , int);
/*----------------------------------------------------------------------------------------------------------------------------------------- */
/*----------------------------------------------------------------------------------------------------------------------------------------- */
/*----------------------------------------------------------------------------------------------------------------------------------------- */
void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{  
	unsigned char *im;
	const int *dimim;
	double *I;
	double *feat , *des;
	int ny , nx , nynx , dimcolor = 1 , npts , i , d;
	double	scale_default[1]    = {1};
	struct opts options;
	mxArray *mxtemp;
	double *tmp , temp , scalemin = 1.0;
	int tempint , dummy;

	options.nscale        = 1;
	options.sigma_scale   = 0.6;
	options.deltax        = 10;
	options.deltay        = 10;
    options.color         = 0;
	options.nbins         = 32;
	options.patchsize     = 16;
	options.norm          = 1;
	options.clamp         = 0.2;

	/* Input 1  */

	if( (mxGetNumberOfDimensions(prhs[0]) > 3) || !mxIsUint8(prhs[0]) )
	{
		mexErrMsgTxt("I must be (ny x nx x [3]) in UINT8 format\n");
	}


	im          = (unsigned char *)mxGetData(prhs[0]);
	dimim       = mxGetDimensions(prhs[0]);
	ny          = dimim[0];
	nx          = dimim[1];
	nynx        = ny*nx;

	/* Input 2  */

	if ((nrhs > 1) && !mxIsEmpty(prhs[1]) )
	{
		mxtemp                             = mxGetField( prhs[1], 0, "scale" );
		if(mxtemp != NULL)
		{
			if((mxGetM(mxtemp) != 1)  && (mxGetN(mxtemp) != 2) )
			{
				mexErrMsgTxt("options.scale must be (1 x nscale) in double format");
			}

			options.scale                  = mxGetPr(mxtemp);
			options.nscale                 = (int) mxGetN(mxtemp);
			for (i = 0 ; i < options.nscale ; i++)
			{
				if((options.scale[i] < 0.0) || (options.scale[i] > 1.0))
				{
					mexErrMsgTxt("0.0 <= options.scale(i)  <= 1.0 ");
				}
				scalemin = min(scalemin , options.scale[i]);
			}
		}
		else
		{
			options.scale                 = (double *)mxMalloc(1*sizeof(double));
			for(i = 0 ; i < options.nscale ; i++)
			{
				options.scale[i]          = scale_default[i];
			}	
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "sigma_scale");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			temp                          = tmp[0];
			if( (temp < 0.0) )
			{
				mexPrintf("sigma_scale > 0 ");	
				options.sigma_scale       = 0.6;

			}
			else
			{
				options.sigma_scale       = temp;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "deltax");
		if(mxtemp != NULL)
		{	
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			dummy                         = (int) floor(nx*scalemin);
			if( (tempint < 1) || (tempint > dummy) )
			{
				mexPrintf("deltax must be > 0 and floor(nx*min(scale)), force to 1");	
				options.deltax           = dummy;
			}
			else
			{
				options.deltax         = tempint;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "deltay");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			dummy                         = (int) floor(ny*scalemin);

			if( (tempint < 1))
			{
				mexPrintf("deltay must be > 0 and floor(ny*min(scale)), force to 1");	
				options.deltay           = dummy;
			}
			else
			{
				options.deltay         = tempint;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "color");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			if( (tempint < 0) || (tempint > 6))
			{
				mexPrintf("color = {0,1,2,3,4,5}, force to 0\n");	
				options.color             = 0;
			}
			else
			{
				options.color            = tempint;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "norm");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			if( (tempint < 0) || (tempint > 4) )
			{
				mexPrintf("norm = {0 , 1 , 2 , 3 , 4}, force to 0");	
				options.norm              = 0;
			}
			else
			{
				options.norm              = tempint;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "clamp");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			temp                          = tmp[0];

			if( (temp < 0.0) )
			{
				mexPrintf("clamp must be >= 0, force to 0.2");	
				options.clamp             = 0.2;
			}
			else
			{
				options.clamp             = temp;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "nbins");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			if( (tempint < 1))
			{
				mexPrintf("nbins > 0, force to 4");	
				options.nbins            = 32;
			}
			else
			{
				options.nbins            = tempint;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "patchsize");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			if( (tempint < 1))
			{
				mexPrintf("patchsize > 0, force to nbins/2");	
				options.patchsize         = options.nbins/2;
			}
			else
			{
				options.patchsize         = tempint;
			}
		}
	}
	else
	{
		options.scale                 = (double *)mxMalloc(options.nscale*sizeof(double));
		for(i = 0 ; i < options.nscale ; i++)
		{
			options.scale[i]          = scale_default[i];
		}	
	}
	if((mxGetNumberOfDimensions(prhs[0]) == 2))
	{
		options.color  = 0;
		I              = (double *)mxMalloc(nynx*sizeof(double));
		dimcolor       = 1;
		for (i = 0 ; i < nynx ; i++)
		{
			I[i]       = (double)im[i];
		}	
	}
	else
	{
		if((options.color == 0) )
		{
			I        = (double *)mxMalloc(nynx*sizeof(double));
			rgb2gray(im , ny , nx , I);
			dimcolor = 1;
		}
		else if (options.color == 1)
		{
			I        = (double *)mxMalloc(3*nynx*sizeof(double));
			for (i = 0 ; i < 3*nynx ; i++)
			{
				I[i] = (double)im[i];
			}
			dimcolor = 3;
		}
		else if (options.color == 2)
		{
			I        = (double *)mxMalloc(3*nynx*sizeof(double));
			rgb2nrgb(im , ny , nx , I);
			dimcolor = 3;
		}
		else if (options.color == 3)
		{
			I        = (double *)mxMalloc(3*nynx*sizeof(double));
			rgb2opponent(im , ny , nx , I);
			dimcolor = 3;
		}
		else if(options.color == 4)
		{
			I        = (double *)mxMalloc(2*nynx*sizeof(double));
			rgb2nopponent(im , ny , nx , I);
			dimcolor = 2;
		}
		else if(options.color == 5)
		{
			I        = (double *)mxMalloc(nynx*sizeof(double));
			rgb2hue(im , ny , nx , I);
			dimcolor = 1;
		}
	}


	/*----------------------- Outputs -------------------------------*/

	d                  =  options.nbins;
	npts               =  options.deltay*options.deltax*options.nscale*dimcolor;

	plhs[0]            =  mxCreateDoubleMatrix(d , npts , mxREAL);
	feat               =  mxGetPr(plhs[0]);

	plhs[1]            =  mxCreateDoubleMatrix(7 , npts , mxREAL);
	des                =  mxGetPr(plhs[1]);


	/*------------------------ Main Call ----------------------------*/

	denseCOLOR(I , ny , nx , options , feat , des , npts , dimcolor);

	/*--------------------------- Free memory -----------------------*/

	if ( (nrhs > 1) && !mxIsEmpty(prhs[1]) )
	{
		if ( (mxGetField( prhs[1] , 0 , "scale" )) == NULL )
		{
			mxFree(options.scale);
		}
	}
	else
	{
		mxFree(options.scale);
	}
	mxFree(I);
}

/*----------------------------------------------------------------------------------------------------------------------------------------- */
void denseCOLOR(double *I , int ny , int nx , struct opts options , double *feat , double *des  , int npts , int dimcolor)
{
	double *scale = options.scale ;
	int nscale = options.nscale , color = options.color;
	int patchsize = options.patchsize , deltax = options.deltax , deltay = options.deltay;
	int nbins = options.nbins , norm = options.norm , d = nbins ;
	double clamp = options.clamp , sigma_scale = options.sigma_scale;
	int s , i , j , x , y , v , nycurrent , nxcurrent , nymax , nxmax , nynxmax , half_patch , half_patch1, nynx = ny*nx , vnynx;
	int hmax, nmax  , indi , indj , indx  , indf , indd;
	int nynew , nxnew  , lynew , lxnew , offsety , offsetx , co , cof , idx , patchsizecurrent;
	double scalecurrent , scalemax = 0.0  , scalemin = 1.0, sigmamax , tmp , sum , mini , maxi , delta , ratio;
	double *kernel , *I_filtered, *aux;

	for (s = 0 ; s < nscale ; s++)
	{
		scalemax                     = max(scalemax , scale[s]);
		scalemin                     = min(scalemin , scale[s]);
	}

	nymax                            = (int) floor(ny*scalemax);
	nxmax                            = (int) floor(nx*scalemax);
	nynxmax                          = nymax*nxmax;

	sigmamax                         = scalemin < 1.0 ? sigma_scale / scalemin : sigma_scale;
	hmax                             = (int) ceil( sigmamax * CTE );
	nmax                             = 1 + 2*hmax;

	kernel                           = (double *)malloc(nmax*sizeof(double));
	I_filtered                       = (double *)malloc(nynxmax*sizeof(double));
	aux                              = (double *)malloc(nymax*nx*sizeof(double));


	co                               = 0;
	for (v = 0 ; v < dimcolor ; v++)
	{
		vnynx    = v*nynx;
		if( (color == 0) || (color == 1) )
		{
			mini     = 0.0;
			maxi     = 255.0;
		}
		else if(color == 2)
		{
			mini     = 0.0;
			maxi     = 1.0;
		}
		else if(color == 3)
		{
			if(v == 0)
			{
				mini     = -255*invsqrt2;
				maxi     = 255*invsqrt2;
			}
			else if(v==1)
			{
				mini     = -510*invsqrt6;
				maxi     = 510*invsqrt6;
			}
			else
			{
				mini     = 0.0;
				maxi     = 765;
			}
		}
		else if(color == 4)
		{
			if(v == 0)
			{
				mini     = -sqrt3/invsqrt2;
				maxi     = sqrt3/invsqrt2;
			}
			else
			{
				mini     = -2.0/invsqrt2;
				maxi     = 1.0/invsqrt2;
			}
		}
		else
		{
			mini     = -PIhalf2;
			maxi     = PIhalf2;
		}

		delta    = (nbins + tiny)/((maxi-mini));

		for (s = 0 ; s < nscale ; s++)
		{
			scalecurrent                     = scale[s];
			patchsizecurrent                 = max(1,(int)floor(patchsize*scalecurrent));
			half_patch                       = (patchsizecurrent - 1)/2;
			half_patch1                      = -half_patch + 1;
			ratio                            = scalemax/scalecurrent;

			if(scalecurrent != 1.0)
			{
				nycurrent                    = (int) floor(ny*scalecurrent);
				nxcurrent                    = (int) floor(nx*scalecurrent);
				gaussian_sampler(I + vnynx , ny , nx , scalecurrent, sigma_scale , kernel , I_filtered , aux);
			}
			else
			{
				nycurrent                    = ny;
				nxcurrent                    = nx;
				for(i = 0 ; i < nynx ; i++)
				{
					I_filtered[i]            = (double) I[i + vnynx];
				}
			}

			nynew        = max(1 , nycurrent - patchsizecurrent);
			nxnew        = max(1 , nxcurrent - patchsizecurrent);
			lynew        = (int) floor(nynew/(deltay - 1) );
			lxnew        = (int) floor(nxnew/(deltax - 1) );
			offsety      = max(0 , (int) ( floor((nynew - (deltay - 1)*lynew)/2))) + half_patch - 1;
			offsetx      = max(0 , (int) ( floor((nxnew - (deltax - 1)*lxnew)/2))) + half_patch - 1;

			for(i = 0 ; i < deltax ; i++)
			{
				indi = offsetx + i*lxnew ;

				for (j = 0 ; j < deltay ; j++)
				{
					indj  = offsety + j*lynew;
					indf  = co*d;
					indd  = co*7;	
					cof   = indf;

					for(x = 0 ; x < patchsizecurrent ; x++)
					{
						indx = (half_patch1 + x + indi)*nycurrent + indj;
						for(y = 0 ; y < patchsizecurrent ; y++)
						{
							tmp       = I_filtered[y + half_patch1 + indx] - mini;
							idx       = (int) floor(tmp*delta) ;
							feat[idx + indf]++;
						}
					}

					for(x = cof ; x < (cof + d) ; x++)
					{
						feat[x] *= ratio;
					}


					tmp  = 0.0;
					if(norm == 1)
					{
						cof   = indf;
						for(x = cof ; x < (cof + d) ; x++)
						{
							tmp += feat[x];
						}

						sum = 1.0/(tmp + tiny);
						for(x = cof ; x < (cof + d) ; x++)
						{
							feat[x] *= sum;
						}
					}
					if(norm == 2)
					{
						cof   = indf;
						for(x = cof ; x < (cof + d) ; x++)
						{
							tmp   += (feat[x]*feat[x]);
						}

						sum  = 1.0/sqrt(tmp + verytiny);
						for(x = cof ; x < (cof + d) ; x++)
						{
							feat[x]  *= sum;
						}
					}

					if(norm == 3)
					{
						cof   = indf;
						for(x = cof ; x < (cof + d) ; x++)
						{
							tmp   += feat[x];
						}

						sum  = 1.0/(tmp + tiny);
						for(x = cof ; x < (cof + d) ; x++)
						{
							feat[x]  *= sqrt(feat[x]*sum);
						}
					}

					if(norm == 4)
					{
						cof   = indf;
						for(x = cof ; x < (cof + d) ; x++)
						{
							tmp   += (feat[x]*feat[x]);
						}

						sum  = 1.0/sqrt(tmp + verytiny);
						for(x = cof ; x < (cof + d) ; x++)
						{
							feat[x]  *= sum;
							if(feat[x] > clamp)
							{
								feat[x] = clamp;
							}
						}
						sum  = 0.0;
						for(x = cof ; x < (cof + d) ; x++)
						{
							sum   += (feat[x]*feat[x]);
						}
						sum  = 1.0/sqrt(sum + verytiny);
						for(x = cof ; x < (cof + d) ; x++)
						{
							feat[x]  *= sum;
						}
					}

					des[0 + indd] = (double)(indj + 1);
					des[1 + indd] = (double)(indi + 1);
					des[2 + indd] = scalecurrent;
					des[3 + indd] = (double) (v+1);
					des[4 + indd] = (double) (nycurrent);
					des[5 + indd] = (double) (nxcurrent);
					des[6 + indd] = (double) (patchsizecurrent);
					co++;
				}
			}
		}
	}
	free(I_filtered);
	free(kernel);
	free(aux);
}
/*----------------------------------------------------------------------------------------------------------------------------------------- */
/*----------------------------------------------------------------------------------------------------------------------------------------- */
/*----------------------------------------------------------------------------------------------------------------------------------------- */
/*----------------------------------------------------------------------------------------------------------------------------------------- */
void gaussian_sampler( double *I  , int ny , int nx , double scalecurrent, double sigma_scale, double *kernel , double *I_filtered , double *aux)
{
	int h , n;
	int nynew ,  nxnew , indx , indxnew ;
	int x,y,xc,yc,i,j;
	double xx,yy,sum;
	int nx2 = 2 * nx;
	int ny2 = 2 * ny;
	double sigmacurrent , tempscale = 1.0/scalecurrent;

	nynew                        = (int) floor(ny*scalecurrent);
	nxnew                        = (int) floor(nx*scalecurrent);
	sigmacurrent                 = scalecurrent < 1.0 ? sigma_scale / scalecurrent : sigma_scale;
	h                            = (int) ceil( sigmacurrent * CTE );
	n                            = 1 + 2*h;

	gaussian_kernel( kernel , sigmacurrent , h , n );

	for(y=0 ; y < nynew ; y++)
	{
		yy      = (double) y * tempscale;          
		yc      = (int) floor( yy + 0.5 ); 

		indx    = 0;
		indxnew = 0;

		for(x=0 ; x < nx ;x++)
		{
			sum      = 0.0;
			for(i = 0 ; i < n ; i++)
			{
				j                = yc - h + i;
				while(j<0)    j +=ny2;
				while(j>=ny2) j -=ny2;
				if(j>=ny)     j  = ny2-1-j;
				sum             += I[ j + indx ] * kernel[i];

			}
			aux[y + indxnew ] = sum;
			indx             += ny;
			indxnew          += nynew;
		}
	}

	indx =  0;
	for(x=0 ; x < nxnew ; x++)
	{
		xx   = (double) x * tempscale;           
		xc   = (int) floor( xx + 0.5 ); 
		for(y=0 ; y < nynew ; y++)
		{
			sum = 0.0;
			for(i=0 ; i < n ; i++)
			{
				j                = xc - h + i;
				while(j<0)    j +=nx2;
				while(j>=nx2) j -=nx2;
				if(j>=nx)     j  = nx2-1-j;
				sum             += aux[ y + j * nynew ] * kernel[i];
			}
			I_filtered[ y + indx ] = sum;
		}
		indx    += nynew;
	}
}

/*----------------------------------------------------------------------------------------------------------------------------------------- */
/*----------------------------------------------------------------------------------------------------------------------------------------- */
/*----------------------------------------------------------------------------------------------------------------------------------------- */
void gaussian_kernel( double *kernel, double sigma , int h , int n )
{
	int i;
	double val , temp_sigma = 1.0/sigma;
	double sum = 0.0;

	/* compute gaussian kernel */
	for(i = 0 ; i < n ; i++)
	{
		val       = (double)( i - h ) * temp_sigma;
		kernel[i] = exp( -0.5 * val * val );
		sum      += kernel[i];
	}

	sum = 1.0/sum;
	/* normalization */
	for(i=0 ; i < n ; i++) 
	{
		kernel[i] *= sum;
	}
}
/*----------------------------------------------------------------------------------------------------------------------------------------- */
/*----------------------------------------------------------------------------------------------------------------------------------------- */
/*----------------------------------------------------------------------------------------------------------------------------------------- */
void rgb2gray(unsigned char *rgb , int Ny , int Nx , double *gray)
{
	/* T    = inv([1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703]);

	vect = T(1 , :);

	*/
	int i , NyNx = Ny*Nx , NyNx2 = 2*NyNx;

	for (i = 0 ; i  < NyNx ; i++)
	{
		gray[i]   = 0.298936021293776*rgb[i] + 0.587043074451121*rgb[i + NyNx] + 0.114020904255103*rgb[i + NyNx2];
	}
}
/*--------------------------------------------------------------------------------------------------------------------------------------------*/
void rgb2nrgb(unsigned char *rgb , int Ny , int Nx , double *nrgb)
{
	int i , NyNx = Ny*Nx , NyNx2 = 2*NyNx ;
	double r , g , b , cte;

	for (i = 0 ; i  < NyNx ; i++)
	{
		r                = (double)rgb[i];
		g                = (double)rgb[i + NyNx];
		b                = (double)rgb[i + NyNx2];
		cte              = (r + g + b);
		if(cte != 0.0)
		{
			cte          = 1.0/cte;
		}
		else
		{
			cte          = 1.0;
		}

		nrgb[i]          = r*cte;
		nrgb[i + NyNx]   = g*cte;
		nrgb[i + NyNx2]  = b*cte;
	}
}
/*--------------------------------------------------------------------------------------------------------------------------------------------*/
void rgb2opponent(unsigned char *rgb , int Ny , int Nx , double *opponent)
{
	int i , NyNx = Ny*Nx , NyNx2 = 2*NyNx ;
	double r , g , b ;

	for (i = 0 ; i  < NyNx ; i++)
	{
		r                    = (double)rgb[i];
		g                    = (double)rgb[i + NyNx];
		b                    = (double)rgb[i + NyNx2];

		opponent[i]          = (r - g)*invsqrt2;
		opponent[i + NyNx]   = (r + g - 2*b)*invsqrt6;
		opponent[i + NyNx2]  = (r + g + b)*invsqrt3;
	}
}
/*--------------------------------------------------------------------------------------------------------------------------------------------*/
void rgb2nopponent(unsigned char *rgb , int Ny , int Nx , double *nopponent)
{
	int i , NyNx = Ny*Nx , NyNx2 = 2*NyNx;
	double r , g , b , cte ;

	for (i = 0 ; i  < NyNx ; i++)
	{
		r                     = (double)rgb[i];
		g                     = (double)rgb[i + NyNx];
		b                     = (double)rgb[i + NyNx2];
		cte                   = (r + g + b);
		if(cte != 0.0)
		{
			cte               = 1.0/(cte*invsqrt2);
		}
		else
		{
			cte               = 1.0;
		}
		nopponent[i]          = sqrt3*(r - g)*cte;
		nopponent[i + NyNx]   = (r + g - 2*b)*cte;
	}
}
/*--------------------------------------------------------------------------------------------------------------------------------------------*/
void rgb2hue(unsigned char *rgb , int Ny , int Nx , double *hue)
{
	int i , NyNx = Ny*Nx , NyNx2 = 2*NyNx ;
	double r , g , b ;

	for (i = 0 ; i  < NyNx ; i++)
	{
		r         = (double)rgb[i];
		g         = (double)rgb[i + NyNx];
		b         = (double)rgb[i + NyNx2];
		hue[i]    = atan(sqrt3*(r - g)/(r + g - 2*b + verytiny));
	}
}
/*--------------------------------------------------------------------------------------------------------------------------------------------*/