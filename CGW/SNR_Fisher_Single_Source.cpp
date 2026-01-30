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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//FUNCTION ANGLES EVOLUTION:
my_type l(my_type t, const T1 &GW_par);  //polarization angle

//ANTENNA PATTERN FUNCTION
my_type F_plus(const T1 &GW_par, const my_type phi_Pi, const my_type theta_Pi);
my_type F_cross(const T1 &GW_par, const my_type phi_Pi, const my_type theta_Pi);

//RESIDUALS:
my_type Sn_plus(my_type t, int n, const T1 &GW_par);
extern "C" void fortran_bessel_jn_(int* n1, int* n2, double* x, double* y);
my_type Sn_cross(my_type t, int n ,const T1 &GW_par);
extern "C" void fortran_bessel_jn_(int* n1, int* n2, double* x, double* y);

//n residual
my_type sn(my_type t, int n ,const T1 &GW_par,const my_type phi_Pi, const my_type theta_Pi);
//Pulsar noise spectral density
my_type Pulsar_PSD(my_type f, const T1 &Pulsar_Noise_Prop);

//peak harmonic emison
my_type n_peak(my_type e);


//SNR function
my_type S_2_fun(my_type t, const T1 &GW_par ,const T1 &Pulsar_Noise_Prop, const my_type phi_Pi,const my_type theta_Pi);

// Define wrapper function for integration
struct S_2_fun_params {

	T1 GW_par;
    T1 Pulsar_Noise_Prop;
    my_type phi_Pi;
    my_type theta_Pi;
};

my_type S_2_fun_wrapper(my_type t, void *params_void) {
    S_2_fun_params *params = (S_2_fun_params *)params_void;
    return S_2_fun(t, params->GW_par, params->Pulsar_Noise_Prop, params->phi_Pi, params->theta_Pi);
};



//funtion to compute SNR over all the pulsar:
my_type compute_SNR(const T1 &GW_par, const T1 &phi_PULSAR_200, const T1 &theta_PULSAR_200, my_type S_white, const T1 &A_rd_200, const T1 &Gamma_rd_200, my_type T_obs_s); 
void printSNR(ofstream& out_file_SNR, const my_type SNR_threshold, T1 &MBHB_properties, const int N_val, int sp);



//Define funtion for interatino succes:
// Custom GSL error handler to prevent program abortion
bool gsl_error_occurred = false; // Flag to track if an error occurred
void my_error_handler(const char *reason, const char *file, int line, int gsl_errno) {
    std::cerr << "GSL Error: " << reason << " (file: " << file << ", line: " << line
              << ", errno: " << gsl_strerror(gsl_errno) << ")" << std::endl;
    gsl_error_occurred = true; // Set the error flag
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
															//START MAIN
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{



	// Catalogue already sorted by Amplitude:
	string path = "/home/rtruant/FedericaPaper/NEW_Catalogs/Catalogs/Cat_GW_Signal/";
	string model_name = argv[1]; //MBHB_Resolved.dat

	//if the file is not foud break:
	ifstream read_model_control(path+model_name);
	if ( !read_model_control.is_open() )
	{
		cout << "Error: file "+path+model_name+" not found!" << endl;
		exit(-1);
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
															//OUTPUT File name
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Amplitude cut
	//my_type A_cut = 6.e-17;
    //SNR treshold
	my_type SNR_threshold = 0.0;

	//SNR outfile path:
	string path_Out = "/home/rtruant/FedericaPaper/NEW_Catalogs/SNR_GWB_SKA/OUTPUT/";
	ofstream out_file_SNR(path_Out+"SNR_"+model_name+".dat");

    // track execution time
	timeval start, end;
	gettimeofday(&start, NULL);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										//Frequency domain and computation lum distance
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    
	//define frequency bins
	cout << "G = " << G << ", c = " << c << endl;
	my_type min_f = 1./(T_obs*year), max_f = 1e-7; // Hz
	const int N_bin_f = int( ( max_f-min_f ) / ( 1./(T_obs*year) ) ) + 1;

	//frequency grid, gwb, and the gwb with the subtracted contribution if the snr of the source is larger than the setted treshold
	T1 F_obs(N_bin_f+1);

	// frequency range at which gwb will be evaluated
	for (int i = 0; i <= N_bin_f; ++i)
	{
		F_obs[i] = min_f + 1./(T_obs*year)*i;
	}
	my_type step_fobs = F_obs[1]-F_obs[0];

	cout<< min_f << ", max_fobs = " << max_f << ", delta_fobs = " << step_fobs << endl;


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


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
															//Read Pulsars postion and noise properties
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //PULSAR SECTION --> read at the beginning and the selct 200 pulsar randomly
	//mean cadence time two weeks:
	my_type t_cad=60.*60.*24.*14.;
	//only white noies:
	my_type sig = 100.e-9;
	//Puslar white NPSD:
	const my_type S_white = 2.*sig*sig*t_cad;
 
	//READ THE PULSARS SKY LOCATION and noise properties: 
	string path_Pulsar = "/home/rtruant/FedericaPaper/NEW_Catalogs/SKA_PTA/Pulsar_Pop/";
	string model_Pulsar = argv[2];
 
	int N_Pulsar = 0;  //number of pulsar considered
	string line_Pulsar;
	ifstream read_model_Pulsar(path_Pulsar+model_Pulsar);
	if ( !read_model_Pulsar.is_open() )
	{
	    cout << "Error: file "+path_Pulsar+model_Pulsar+" not found!" << endl;
		exit(-1);
	}
	while (!read_model_Pulsar.eof())
	{
		getline(read_model_Pulsar, line_Pulsar);
		N_Pulsar++;
	}
	read_model_Pulsar.close();
 
	N_Pulsar = N_Pulsar - 1; //last line is read two times
	cout<<"The tolal number of pulsar is Np= "<<N_Pulsar<<endl;
 
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
	//open the file and read RA DEC
	read_model_Pulsar.open(path_Pulsar+model_Pulsar);
	read_model_Pulsar.precision(18);
	//noise properties:
	my_type phi_P, theta_P,A_rd_P,Gamma_rd_P;
	T1 phi_PULSAR_200(N_Pulsar), theta_PULSAR_200(N_Pulsar), A_rd_200(N_Pulsar), Gamma_rd_200(N_Pulsar);
	 
	for (int n_p = 0; n_p < N_Pulsar; n_p++)
	{  
		read_model_Pulsar >> phi_P >> theta_P>>A_rd_P>>Gamma_rd_P;
		phi_PULSAR_200[n_p] = phi_P;    //phi [0,2pi] 
		theta_PULSAR_200[n_p] = theta_P; //theta [0,pi] spherical polar coordinates
		A_rd_200[n_p] = A_rd_P;
		Gamma_rd_200[n_p] = Gamma_rd_P;
	}
	read_model_Pulsar.close();
 
	//end pulsar section


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
											//Read MBHB properties
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//READ THE SOURCE PARAMETERS:
	// find the number of sources by reading the lines of the input file and check path
	int N_sources = countLinesInFile(path,model_name);
	
	cout<<"The number of source is N = "<<N_sources<<endl;
	//////////////////////////////////////////////////////////////////////////////
	// open the file to read source parameters


	//binary system properties:
	my_type Mc_5_3;
	my_type M1_loc, M2_loc, z_loc, f_gw_p_loc, e_loc, fk_r_loc, i_loc, PSI_loc, l0_loc, y0_loc, RA_loc, DEC_loc;


	//start write outfile 
	out_file_SNR.precision(10);   
	ifstream read_model(path+model_name);
	read_model.precision(18);
	  
	//calculate the srn only for the possible resolvable source--->i  want to reset the reading only after  res source
	for(int k=0; k < N_sources; k++)
	{
	  
		read_model >> M1_loc >> M2_loc >> z_loc >> f_gw_p_loc >> e_loc  >> fk_r_loc  >> i_loc  >> PSI_loc >> l0_loc >> y0_loc >> RA_loc  >> DEC_loc; 

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
												//Defining Quantities to compute SNR and Fisher 
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//rest frame chirp mass
		my_type Mc_5_3 = (M1_loc*M2_loc)/(pow(M1_loc+M2_loc,1./3.));

		//lets calculate the orbital observed frequency self consistently; 
	  	my_type f_r = fk_r_loc; //restframe orb freq (Hz)
		my_type f = f_r / (1.+ z_loc); //obs keplerian freq
		my_type f_2 = f*2.;//obs GW freq second harmonic

		//variables i need to calculate the GWB and the noise spectral density
	    my_type g_n_e, f_n, Si_n, Sh_n; 
                 
		if (e_loc ==0) e_loc=0.00001; //stability fisher matrix derivative
		//harmonic at which the emission is higher
		my_type n_max = 4.*n_peak(e_loc);

		//Calculate the comoving distance
		int id = int( (z_loc - min_zD)/delta_zD );
		my_type xd = (z_loc-zD_array[id]) / (zD_array[id+1]-zD_array[id]);
		my_type D_M_distance = comoving_distance_array[id]*(1-xd) + comoving_distance_array[id+1]*xd; //Mpc

		//from deg we want to go in radians:
		my_type phi_RAD = RA_loc; // phi [0,2pi] already in rad
	    my_type theta_RAD = DEC_loc;  //theta [0,pi] spherical polar coordinates

		//all quantities are observed:
		my_type Mc_5_3_z = Mc_5_3*pow(1.+z_loc,5./3.); //redshifted chirp mass
        my_type DL = D_M_distance * (1.+z_loc);        //luminosity distance
        my_type C = Mc_5_3_z / DL * G_5_3 / c_4;       //Residual amplitude

		T1 GW_par = { C , f, e_loc, i_loc, PSI_loc, l0_loc, y0_loc, phi_RAD, theta_RAD}; //9 free parameters


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
														//Compute SNR and print results 
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//COMPUTE THE SNR:
		my_type SNR = compute_SNR(GW_par,phi_PULSAR_200, theta_PULSAR_200, S_white ,A_rd_200, Gamma_rd_200, T_obs_s);
				
		cout<<"SNR first = "<<SNR<<endl;

		//PRINT S/N AND SOURCE PROPERTIES:
		int N_val = 12; //Number of quantites to print
		//Array contating the properties to print
		T1 MBHB_Properties = { SNR, M1_loc, M2_loc, f_r, e_loc, z_loc, i_loc, PSI_loc, RA_loc, DEC_loc, l0_loc, y0_loc}; //12
		printSNR(out_file_SNR, SNR_threshold, MBHB_Properties, N_val, sp ); 

	}

	//close outfiles:
	out_file_SNR.close();
	read_model.close();
		  
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


//FUNCTION ANGLES EVOLUTION:

//Orbital phase
my_type l(my_type t, const T1 &GW_par){
	my_type l0 = GW_par[5];   //intial phase
	my_type f = GW_par[1];    //observed keplerian freq

	return l0 + 2.*pi*f*t;
}

//ANTENNA PATTERN FUNCTIONS:
my_type F_plus(const T1 &GW_par, const my_type phi_Pi, const my_type theta_Pi){
	//Rosado et al 2015, k and j polarization basis vector 
	my_type phi = GW_par[7];   //right ascension GW source RAD
	my_type theta = GW_par[8];  //declination GW sourve RAD
	my_type PSI_GW =GW_par[4];         //Polarization angle

	T1 u = { cos(PSI_GW)*cos(theta)*cos(phi) - sin(PSI_GW)*sin(phi),  //polarization basis tensor
             cos(PSI_GW)*cos(theta)*sin(phi) + sin(PSI_GW)*cos(phi),
            -cos(PSI_GW)*sin(theta)}; 

    T1 v = { sin(PSI_GW)*cos(theta)*cos(phi) + cos(PSI_GW)*sin(phi),  //polatization basis tensor
             sin(PSI_GW)*cos(theta)*sin(phi) - cos(PSI_GW)*cos(phi),
            -sin(PSI_GW)*sin(theta)}; 

    T1 Om = {(-sin(theta)*cos(phi)),(-sin(theta)*sin(phi)),(-cos(theta))}; //GW propagation direction in polar spherical coord

    T1 p = {sin(theta_Pi)*cos(phi_Pi),sin(theta_Pi)*sin(phi_Pi),cos(theta_Pi)}; //Pulsar sky location in polar spherical coord
    

    my_type u_p =  u[0]*p[0]+u[1]*p[1]+u[2]*p[2];
    my_type v_p =  v[0]*p[0]+v[1]*p[1]+v[2]*p[2];
    my_type Om_p = Om[0]*p[0]+Om[1]*p[1]+Om[2]*p[2];

    return (1./2.)*((u_p*u_p)-(v_p*v_p))/(1.+Om_p);
}

my_type F_cross(const T1 &GW_par,const  my_type phi_Pi,const my_type theta_Pi ) 
{
    //Rosado et al 2015, k and j polarization basis vector 
	my_type phi = GW_par[7];   //right ascension GW source RAD
	my_type theta = GW_par[8];  //declination GW sourve RAD
	my_type PSI_GW =GW_par[4];         //Polarization angle

	T1 u = { cos(PSI_GW)*cos(theta)*cos(phi) - sin(PSI_GW)*sin(phi),  //polarization basis tensor
             cos(PSI_GW)*cos(theta)*sin(phi) + sin(PSI_GW)*cos(phi),
            -cos(PSI_GW)*sin(theta)}; 

    T1 v = { sin(PSI_GW)*cos(theta)*cos(phi) + cos(PSI_GW)*sin(phi),  //polatization basis tensor
             sin(PSI_GW)*cos(theta)*sin(phi) - cos(PSI_GW)*cos(phi),
            -sin(PSI_GW)*sin(theta)}; 

    T1 Om = {(-sin(theta)*cos(phi)),(-sin(theta)*sin(phi)),(-cos(theta))}; //GW propagation direction in polar spherical coord

    T1 p = {sin(theta_Pi)*cos(phi_Pi),sin(theta_Pi)*sin(phi_Pi),cos(theta_Pi)}; //Pulsar sky location in polar spherical coord
    

    my_type u_p =  u[0]*p[0]+u[1]*p[1]+u[2]*p[2];
    my_type v_p =  v[0]*p[0]+v[1]*p[1]+v[2]*p[2];
    my_type Om_p = Om[0]*p[0]+Om[1]*p[1]+Om[2]*p[2];

    return ((u_p)*(v_p))/(1.+Om_p);
}


//RESIDUALS

my_type Sn_plus(my_type t, int n ,const T1 &GW_par){
	int n1 = n-2, n2= n+2;
    double n_e = n*GW_par[2];
    double jn[5];
    fortran_bessel_jn_(&n1,&n2,&n_e,jn);//between n-2,n+2, argument ne
	// cout << jn[0] << " " << jn[1] << " " << jn[2] << " " << jn[3] << " " << jn[4] << endl;
    //         n-2,            n-1             n                n+1            n+2

    //Prameters:
	my_type C  = GW_par[0]; 			//residual amplitude
	my_type f  = GW_par[1];           //observed keplerian frequency
	my_type w  = 2.*M_PI*f;         //observed keplerian angular velocity
	my_type e  = GW_par[2];           //eccentricity
    my_type i  = GW_par[3];           //inclination angle  
	my_type y0 = GW_par[6];           //pericenter
	
    //functions:
	my_type an =  -C*pow(w,-1./3.) * ( jn[0] - 2.*e*jn[1] + (2./n)*jn[2] + 2.*e*jn[3] - jn[4] ) * sin(n*l(t,GW_par));
	my_type bn =   C*pow(w,-1./3.) * sqrt(1.-e*e)*( jn[0] - 2.*jn[2] + jn[4] ) * cos(n*l(t,GW_par));
	my_type cn =  (2./n)*C*pow(w,-1./3.)*jn[2]*sin(n*l(t,GW_par));

	//function at initial conditions:
	my_type an0 = -C*pow(w,-1./3.) * ( jn[0] - 2.*e*jn[1] + (2./n)*jn[2] + 2.*e*jn[3] - jn[4] ) * sin(n*l(0,GW_par));
	my_type bn0 =  C*pow(w,-1./3.) * sqrt(1.-e*e)*( jn[0] - 2.*jn[2] + jn[4] ) * cos(n*l(0,GW_par));
	my_type cn0 = (2./n)*C*pow(w,-1./3.)*jn[2]*sin(n*l(0,GW_par));

	return (  (-1.)*(1.+cos(i)*cos(i)) * ( an*cos(2.*y0) - bn*sin(2.*y0) ) + (1.-cos(i)*cos(i))*cn   ) - ( ((-1.)*(1.+cos(i)*cos(i)) * ( an0*cos(2.*y0) - bn0*sin(2.*y0) )) + ((1.-cos(i)*cos(i))*cn0) );
}

my_type Sn_cross(my_type t, int n ,const T1 &GW_par){
	int n1 = n-2, n2= n+2;
    double n_e = n*GW_par[2];
    double jn[5];
    fortran_bessel_jn_(&n1,&n2,&n_e,jn);//between n-2,n+2, argument ne
	// cout << jn[0] << " " << jn[1] << " " << jn[2] << " " << jn[3] << " " << jn[4] << endl;
    //         n-2,            n-1             n                n+1            n+2

    //Prameters:
	my_type C  = GW_par[0]; 			//residual amplitude
	my_type f  = GW_par[1];           //observed keplerian frequency
	my_type w  = 2.*M_PI*f;         //observed keplerian angular velocity
	my_type e  = GW_par[2];           //eccentricity
    my_type i  = GW_par[3];           //inclination angle  
	my_type y0 = GW_par[6];           //pericenter
	
    //functions:
	my_type an =  -C*pow(w,-1./3.) * ( jn[0] - 2.*e*jn[1] + (2./n)*jn[2] + 2.*e*jn[3] - jn[4] ) * sin(n*l(t,GW_par));
	my_type bn =   C*pow(w,-1./3.) * sqrt(1.-e*e)*( jn[0] - 2.*jn[2] + jn[4] ) * cos(n*l(t,GW_par));

	//function at initial conditions:
	my_type an0 = -C*pow(w,-1./3.) * ( jn[0] - 2.*e*jn[1] + (2./n)*jn[2] + 2.*e*jn[3] - jn[4] ) * sin(n*l(0,GW_par));
	my_type bn0 =  C*pow(w,-1./3.) * sqrt(1.-e*e)*( jn[0] - 2.*jn[2] + jn[4] ) * cos(n*l(0,GW_par));

	return ( 2.*cos(i)*( bn*cos(2.*y0) + an*sin(2.*y0) ))  -  ( 2.*cos(i)*( bn0*cos(2.*y0) + an0*sin(2.*y0) ) );
}


//n residuals:
my_type sn(my_type t, int n ,const T1 &GW_par,const my_type phi_Pi,const my_type theta_Pi){
	return ((F_plus(GW_par,phi_Pi,theta_Pi)*Sn_plus(t,n,GW_par)) + 
		(F_cross(GW_par,phi_Pi,theta_Pi)*Sn_cross(t,n,GW_par))) ;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//n_peak
my_type n_peak(my_type e){

	my_type c1 = -1.01678, c2 = 5.57372, c3 = -4.9271, c4 = 1.68506;
	return   2*( 1 + c1*pow(e,1.) + c2*pow(e,2.) +c3*pow(e,3.) +c4*pow(e,4.) ) * pow((1-e*e),-3/2);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//pulsar NPSD:
my_type Pulsar_PSD(my_type f, const T1 &Pulsar_Noise_Prop){

    //Pulsar NPSD
	my_type Sp;

	//Assume for SGWB simple power law EPTA resutls
	my_type Agw = 2.5e-15;
	my_type alpha_gwb = -2./3.;
	my_type y_gwb = 13./3.;

	my_type Sh_gwb = Agw*Agw / (12.*pi*pi*pow(f_year, 2.*alpha_gwb )) * pow(f,-1.*y_gwb);

	//Pulsar noise prop:
	my_type S_w = Pulsar_Noise_Prop[0];
	my_type A_red = Pulsar_Noise_Prop[1];
	my_type y_red = Pulsar_Noise_Prop[2];
	
	//dependig of the fiiting model for the noise PSD compute differently:

	my_type Sh_red = A_red*A_red / (12.*pi*pi) * pow( f / f_year, -1.*y_red )*year*year*year;

	Sp =  S_w + Sh_red + Sh_gwb;
	
	return Sp;
}


//Function to integrate

//S_tot integran dunction
my_type S_2_fun(my_type t, const T1 &GW_par, const T1 &Pulsar_Noise_Prop, const my_type phi_Pi, const my_type theta_Pi){

     

    //GW parameters
	my_type e = GW_par[2];
	my_type f = GW_par[1];
	my_type min_f = 1./T_obs_s;
	my_type step_fobs = 1./T_obs_s;

	//select the frequency bin


	int n_max = 4*ceil(n_peak(e));

	my_type S2=0.;

	for(int n = 1; n <= n_max; n++)
	{   
        //minum freq and number of bins --> change this
		my_type min_f = 1./(T_obs*year), max_f = 1e-7; // Hz
	    const int N_bin_f = int( ( max_f-min_f ) / ( 1./(T_obs*year) ) ) + 1;

		my_type f_n = n*f;    //oberved frame harmonic frequence

		my_type index = (f_n-min_f)/step_fobs;//frequency bin:

        //if signal outside the pta band no signal -->rediual equal zero
		if (index > N_bin_f)   break;
	                    
		//if signal outside the pta band no signal -->rediual equal zero
		if (index < 0) continue;

		//PUlsar NPSD (pulsar + gwb)
		my_type Sp_n = Pulsar_PSD(f_n,Pulsar_Noise_Prop);
						        
		S2 += 2.*sn(t,n ,GW_par, phi_Pi, theta_Pi)*sn(t,n,GW_par, phi_Pi, theta_Pi) / Sp_n;
	}

	
	return S2;
}


//SNR computaion over all the pulsar
my_type compute_SNR(const T1 &GW_par, const T1 &phi_PULSAR_200, const T1 &theta_PULSAR_200, my_type S_white, const T1 &A_rd_200, const T1 &Gamma_rd_200, my_type T_obs_s) {


	my_type SNR2_tot = 0.;
	my_type SNR2_i = 0.;


	// START FOR LOOP ON THE PULSARS:
	for (int n_p = 0; n_p < 200; n_p++) 
	{
		// Reset the SNR_i:
		SNR2_i = 0.;

		// Sky location of the pulsar, already in radians
		my_type phi_Pi = phi_PULSAR_200[n_p];
		my_type theta_Pi = theta_PULSAR_200[n_p];

		// Array with the pulsar properties:
		//               0        1              2                    
		T1 Pulsar_Noise_Prop = {S_white, A_rd_200[n_p], Gamma_rd_200[n_p]};

		// Set up parameters struct
		S_2_fun_params params = {GW_par, Pulsar_Noise_Prop, phi_Pi, theta_Pi};

		// Compute the SNR
		// Set up GSL integration
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);

		// Set up GSL function
		gsl_function F;
		F.function = &S_2_fun_wrapper;
		F.params = &params;

		// Integrate using QAGS routine
		my_type result, error;

		// Set custom error handler
		gsl_set_error_handler(&my_error_handler);

		// Initial tolerances
		my_type epsabs = 1e-2;
		my_type epsrel = 1e-2;

		// Flag to track if we should skip this iteration
		bool skip_iteration = false;

		int status = gsl_integration_qag(&F, 0, T_obs_s, epsabs, epsrel, 1000, GSL_INTEG_GAUSS61, w, &result, &error);
		
		// First attempt at integration
		gsl_error_occurred = false; // Reset error flag

		if (status != GSL_SUCCESS && gsl_error_occurred) 
		{
			// Relax tolerances (e.g., increase by a factor of 10)
			epsabs = 1e-1; // Relaxed absolute tolerance
			epsrel = 1e-1; // Relaxed relative tolerance

			// Retry integration
			gsl_error_occurred = false; // Reset error flag
			status = gsl_integration_qag(&F, 0, T_obs_s, epsabs, epsrel, 1000, GSL_INTEG_GAUSS61, w, &result, &error);

			if (status != GSL_SUCCESS && gsl_error_occurred) 
			{
				skip_iteration = true; // Mark this iteration to be skipped
			}
		}

		// If integration succeeded, update SNR2_i
		if (!skip_iteration)
		{
			SNR2_i = result;
		} 
		else 
		{
			SNR2_i = 0.;
		}

		// Add pulsar SNR to total SNR
		SNR2_tot += SNR2_i;

		// Free GSL workspace
		gsl_integration_workspace_free(w);
		
	} // end loop for on the pulsar


	return sqrt(SNR2_tot); //square root of the SNR

	
}

//Print SNR function:
void printSNR(ofstream& out_file_SNR, const my_type SNR_threshold, T1 &MBHB_properties, const int N_val, int sp) {

    my_type SNR = MBHB_properties[0];
	// Print parameters only if SNR exceeds threshold
	if (SNR > SNR_threshold) 
	{
		for (int i = 0; i < N_val; ++i) 
		{
			out_file_SNR << left << scientific
			<< setw(sp) << MBHB_properties[i];
		}	

		out_file_SNR << endl;
	}

}
