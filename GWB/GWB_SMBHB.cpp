#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <limits>
#include <sys/time.h>
#include <ctime>
#include <chrono>
#include <complex>
#include <omp.h>
#include <random>
#include <algorithm>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include "read_fun.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace std;

typedef double my_type;
typedef long double my_type_bin;

typedef vector< vector< vector< vector<my_type> > > > T4;//tensor 4 indices
typedef vector< vector< vector<my_type> > > T3;//tensor 3 indices
typedef vector< vector<my_type> > T2;//tensor 2 indices -> matrix
typedef vector<my_type> T1;//tensor 1 indices -> vector

typedef vector<int> vec_of_int;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTANTS AND PARAMETERS

int sp = 30;

const my_type Gmks = 6.67384e-11;// m^3/(kg s^2)
const my_type cmks = 299792458.;// m/s
const my_type parsec = 3.08567758e16;// m
const my_type year = 31557600.;// s
const my_type solar_mass = 1.98855e30;// kg
const my_type Mpc = parsec*1e6; // m

// const my_type Mpc_3 = (Mpc*Mpc*Mpc); // m^3
const my_type H0 = 67.3*1000/Mpc;// s^-1
// const my_type mtsun = Gmks*solar_mass/(cmks*cmks*cmks);
// const my_type Rsch = Gmks*solar_mass/(cmks*cmks);

const my_type G = Gmks/pow(Mpc,3) * solar_mass; // Mpc^3/(s^2 Msun)
const my_type c = cmks/Mpc; // Mpc/s
const my_type G_5_3 = pow(G,5./3.), G_10_3 = G_5_3*G_5_3, c_5 = pow(c,5), c_3 = pow(c,3), c_8 = c_5*c_3, c_4=pow(c,4); 
const my_type G_2_3 = pow(G,2./3.), G_4_3 = pow(G,4./3.), G_1_3 = pow(G,1./3.);
	
const my_type pi = M_PI;
const my_type two_pi_2_3 = pow(2*pi,2./3.), two_pi_4_3 = pow(2*pi,4./3.), two_pi_8_3 = pow(2*pi,8./3.), two_pi_10_3 = pow(2*pi,10./3.);

// cosmological parameters
const T1 cosmo_par = {H0, 0., 0.315, 0., 0.683};

// this is the timespan of PTA observations
const my_type f_year = 1./(60.*60.*24.*365.25);
const my_type T_obs =30.; // years
const my_type T_obs_s =T_obs * year;   
const my_type t_cad = 60.*60.*24.*14.;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//FUNCTION
my_type E(my_type z, const T1 &args);
my_type comoving_distance(my_type z, const T1  &args);

my_type g_func(int n, double e);
extern "C" void fortran_bessel_jn_(int* n1, int* n2, double* x, double* y);

//n_peak
my_type n_peak(my_type e);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Catalogue already sorted by Amplitude:
	string path_Models = "SMBHB_populations/";

	//model name + realization name             
	string SMBHB_popName = argv[1];
	    

    // track execution time
	timeval start, end;
	gettimeofday(&start, NULL);
    
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //GWB outfile path
	string out_path_GWB = "OUTPUT/";
	ofstream out_file_GWB(out_path_GWB+"GWB_"+SMBHB_popName);

    // frequency bins and array GWB
    cout << "G = " << G << ", c = " << c << endl;
    my_type min_f = 1./(T_obs*year), max_f = 1e-7; // Hz
    const int N_bin_f = int( ( max_f-min_f ) / ( 1./(T_obs*year) ) ) + 1;

   //frequency grid, gwb, and the gwb with the subtracted contribution if the snr of the source is larger than the setted treshold
    T1 F_obs(N_bin_f+1), hc_gwb_2_vec(N_bin_f+1), hc_gwb_2_vec_sub(N_bin_f+1) ;
	//array containing the SNR and strain for each bin
	T1 SNR_2_bin(N_bin_f+1), hc_2_bin(N_bin_f+1);
  
	// frequency range at which gwb will be evaluated
	for (int i = 0; i <= N_bin_f; ++i)
	{
	    F_obs[i] = min_f + 1./(T_obs*year)*i;
	    //cout << F_obs[i] << endl;
	    hc_gwb_2_vec[i] = 0.;
		hc_gwb_2_vec_sub[i]=0.;
	}
	my_type step_fobs = F_obs[1]-F_obs[0];
    
	cout<< min_f << ", max_fobs = " << max_f << ", delta_fobs = " << step_fobs << endl;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// evaluate array of comoving distances such that during binary sampling the comoving distance can be obtained by interpolation
	int N_zD = 10000;
	my_type min_zD = 0.001, max_zD = 2.5, delta_zD;
	T1 zD_array(N_zD+1), comoving_distance_array(N_zD+1);
	for (int i = 0; i <= N_zD; i++)
	{
	    zD_array[i] = min_zD + (max_zD-min_zD)*i/my_type(N_zD);
	    comoving_distance_array[i] = c/H0 * comoving_distance(zD_array[i], cosmo_par);
	}
	delta_zD = zD_array[2]-zD_array[1];

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//READ THE SOURCE PARAMETERS:
	// find the number of sources by reading the lines of the input file
	int N_sources = 0;
	string line;
	ifstream read_model(path_Models+SMBHB_popName);
	if ( !read_model.is_open() )
	{
	    cout << "Error: file "+path_Models+SMBHB_popName+" not found!" << endl;
	    exit(-1);
	}
	while (!read_model.eof())
	{
	    getline(read_model, line);
	    N_sources++;
	}
	read_model.close();
	//the last source is readed two times:
	N_sources = N_sources - 1;
	cout<<"The number of source is N = "<<N_sources<<endl;
	//////////////////////////////////////////////////////////////////////////////

	// open the file to read source parameters
	read_model.open(path_Models+SMBHB_popName);
	read_model.precision(18);
	if ( !read_model.is_open() )
    {
	    cout << "Error: file "+path_Models+SMBHB_popName+" not found!" << endl;
	    exit(-1);
	}

	//binary system properties:
	my_type Mc_5_3;
	my_type M1_loc, M2_loc, f_gw_obs_loc, a_loc, e_loc, z_loc, i_loc, PSI_loc, RA_loc, DEC_loc, l0_loc, y0_loc;
    //pol contribution
	my_type a,b;
	my_type MeanAng;   

	//vector containing the index of the possible resolvable source
	vec_of_int Id;
	//number of possible resolvable sources:
	int N_res_source=0; 
    
    //calculate the GWB and save the index with amplitude larger than1e-16:
	for (int n_s = 0; n_s < N_sources; n_s++)
	{
	    read_model >> M1_loc >> M2_loc >> f_gw_obs_loc >> a_loc >> e_loc >> z_loc >> i_loc >> PSI_loc >> RA_loc >> DEC_loc >> l0_loc >> y0_loc;
                        
	    /////////////////////////////////////////////////////////////////////
	    //e_loc=0.;
	    //chirp mass
	    Mc_5_3 = pow(10., log10_Mc * 5./3. );
        if (e_loc ==0) e_loc=0.000001; //stability fisher matrix derivative
	    //lets calculate the orbital observed frequency self consistently; 
	    my_type f_kep_rest = (1/(2*math.pi)) * sqrt(  G*(M1_loc+M2_loc)/(a_loc*a_loc*a_loc)); //restframe kep freq (Hz)
	    my_type f_kep_obs = f_kep_rest/(1.+z_loc); //restframe orb freq (Hz)
		my_type f_gw_obs = 2.*f_kep_obs; //gw observed freq (Hz)

        //Derive the comoving distance
		int id = int( (z_loc - min_zD)/delta_zD );
		my_type xd = (z_loc-zD_array[id]) / (zD_array[id+1]-zD_array[id]);
		my_type D_M_distance = comoving_distance_array[id]*(1-xd) + comoving_distance_array[id+1]*xd; //Mpc

		//pol contribution:
		a = 1.+cos(i_loc)*cos(i_loc);
	 	b = -2.*cos(i_loc);
		MeanAng = sqrt( 1./2.* (a*a+b*b) ); 

        //bessel function combination, harmonic freq, gwb values
		my_type g_n_e, f_n, h2_n_square;

        //harmonic at which the emmision is brightest
		my_type n_max = 4*n_peak(e_loc);

		//calculate gw amplitude:
		my_type A_max = 2.*G_5_3*Mc_5_3*(pow(2.*f_gw_obs*pi*(1+z_loc) ,2./3.))/(D_M_distance*c_4); //amplitude; 
    
		// compute the GWB by adding sources
		for (int n = 1; n <= ceil(n_max); n++)// cycle over harmonics
	   	{
		    g_n_e = g_func(n, e_loc);
		    if(g_n_e==0) continue;

		    f_n = n*f_kep_rest;

		    my_type index = (f_n/(1.+z_loc)-min_f)/step_fobs;

		    if (index > N_bin_f) break;
			
			if (index < 0) continue;
			
			// Amaro-seoane et al. 2010
			h2_n_square = 2.*2.*MeanAng*MeanAng* 4.*G_5_3*G_5_3*Mc_5_3*Mc_5_3/(c_8 * D_M_distance*D_M_distance)*pow(2.*pi*f_kep_rest,4./3.) / (n*n) * g_n_e * f_n/(1.+z_loc) / step_fobs;
			hc_gwb_2_vec[int(index)] += h2_n_square; //gwb
			hc_gwb_2_vec_sub[int(index)] += h2_n_square; //gwb

		}

	}//end gwb
    
	read_model.close();
	out_file_GWB.precision(10);

	for (int w = 0; w < N_bin_f; w++)
	{
	    out_file_GWB << left << scientific
	    << setw(sp) << log10( 0.5*(F_obs[w]+F_obs[w+1]))
	    << setw(sp) << log10( sqrt(hc_gwb_2_vec[w]) )
	    << endl;
	}
	out_file_GWB.close();

	cout << " done." << endl;


	gettimeofday(&end, NULL);
	double time_sampling = (double) ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	//save the subtracted gwb:

	cout << "Time elapsed: " << time_sampling << "sec" << endl;	


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	return 0;    	

}






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS

// function entering the comoving distance integral
my_type E(my_type z, const T1 &args)
{
	my_type Omega_r = args[1];
	my_type Omega_M = args[2];
	my_type Omega_k = args[3];
	my_type Omega_L = args[4];

	return 1./sqrt( Omega_r*(1.+z)*(1.+z)*(1.+z)*(1.+z)
								+ Omega_M*(1.+z)*(1.+z)*(1.+z)
								+ Omega_k*(1.+z)*(1.+z)
								+ Omega_L);
}

// function that computes (thorugh trapezoidal rule) the comoving distance up to a redshift z, given the cosmology
my_type comoving_distance(my_type z, const T1 &args)
{
	int NN = 500;
	my_type zz_l, zz_u, D_M = 0;
	for (int l = 0; l < NN; ++l)
	{
		zz_l = z/NN*l;
		zz_u = z/NN*(l+1);
		D_M += ( E(zz_l,args)+E(zz_u,args) ) * (zz_u-zz_l)/2.;
	}

	return D_M;
}

my_type g_func(int n, double e)
{
	int n1 = n-2, n2 = n+2;
	double n_e = n*e;
	double jn[5];
	fortran_bessel_jn_(&n1,&n2,&n_e,jn);

	// cout << jn[0] << " " << jn[1] << " " << jn[2] << " " << jn[3] << " " << jn[4] << endl;

	my_type term1 = (jn[0] - 2.*e*jn[1] + 2./n*jn[2] + 2.*e*jn[3] - jn[4]);
	my_type term2 = (jn[0] - 2.*jn[2] + jn[4]);
	my_type term3 = 4./(3.*n*n)*jn[2]*jn[2];
	my_type n4 = my_type(n);
	n4 *= n4*n4*n4;

	return n4/32.*( term1*term1 + (1.-e*e)*term2*term2 + term3 );
}


//n_peak
my_type n_peak(my_type e){

	my_type c1 = -1.01678, c2 = 5.57372, c3 = -4.9271, c4 = 1.68506;
	return   2*( 1 + c1*pow(e,1.) + c2*pow(e,2.) +c3*pow(e,3.) +c4*pow(e,4.) ) * pow((1-e*e),-3/2);
}