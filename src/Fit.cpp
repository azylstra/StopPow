#include "Fit.h"

struct data {
  size_t n;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> sigma;
};

// ----------------------------------------------------
//				Direct Gaussian Fitting
// ----------------------------------------------------
// Gaussian function, in form to be used with gsl multifit library
int gauss_f(const gsl_vector * p, void * data, gsl_vector * f) {
    // retrieve data:
    size_t n = ((struct data *)data)->n;
    std::vector<double> x = ((struct data *)data)->x;
    std::vector<double> y = ((struct data *)data)->y;
    std::vector<double> stdev = ((struct data *) data)->sigma;

    // fitting parameters:
    double A = gsl_vector_get (p, 0);
    double mu = gsl_vector_get (p, 1);
    double sigma = gsl_vector_get (p, 2);

    // loop over data, putting (gauss - y[i])/sigma[i] into f
    for(size_t i=0; i < n; i++)
    {
        double g = (A/(sqrt(2*M_PI)*sigma)) * exp (-1.*pow(x[i]-mu,2)/(2.*sigma*sigma));
        // sanity check:
        if(stdev[i] > 0)
        	gsl_vector_set(f, i, (g - y[i])/stdev[i] );
        else
        	gsl_vector_set(f, i, 0);
    }

    return GSL_SUCCESS;
}

// Jacobian matrix
int gauss_df(const gsl_vector * p, void * data, gsl_matrix * J) {
    // retrieve data:
    size_t n = ((struct data *)data)->n;
    std::vector<double> x = ((struct data *)data)->x;
    std::vector<double> y = ((struct data *)data)->y;
    std::vector<double> stdev = ((struct data *) data)->sigma;

    // fitting parameters:
    double A = gsl_vector_get (p, 0);
    double mu = gsl_vector_get (p, 1);
    double sigma = gsl_vector_get (p, 2);

    // loop over data, putting (gauss - y[i])/sigma[i] into f
    for(size_t i=0; i < n; i++)
    {
        // Calculate derivatives:
        double dGdA = (1./(sqrt(2*M_PI)*sigma))
            * exp(-1.*pow(x[i]-mu,2)/(2.*sigma*sigma));
        double dGdmu = (A/(sqrt(2*M_PI)*pow(sigma,3)))
            * exp(-1.*pow(x[i]-mu,2)/(2.*sigma*sigma))
            * (x[i] - mu);
        double dGds = (A/(sqrt(2*M_PI)*pow(sigma,2)))
            * exp(-1.*pow(x[i]-mu,2)/(2.*sigma*sigma))
            * ( pow((x[i]-mu)/sigma,2) - 1.0);

        // Set Jacobian, with sanity check:
        if(stdev[i] > 0)
        {
        	gsl_matrix_set (J, i, 0, dGdA/stdev[i]);
        	gsl_matrix_set (J, i, 1, dGdmu/stdev[i]);
        	gsl_matrix_set (J, i, 2, dGds/stdev[i]);
    	}
    	else
    	{
        	gsl_matrix_set (J, i, 0, 0);
        	gsl_matrix_set (J, i, 1, 0);
        	gsl_matrix_set (J, i, 2, 0);
    	}
    }

    return GSL_SUCCESS;
}

// Convenient combination function for above, which GSL needs
int gauss_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
  gauss_f (x, data, f);
  gauss_df (x, data, J);

  return GSL_SUCCESS;
}

// Find the index corresponding to a maximum value in an array
int find_max_i(std::vector<double> & x)
{
    double max = x[0];
    int max_i = 0;
    for(int i=0; i<x.size(); i++)
    {
        if( x[i] > max )
        {
            max = x[i];
            max_i = i;
        }
    }
    return max_i;
}

// Print the state of a multifit
// from GSL example
void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          (unsigned int)iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_blas_dnrm2 (s->f));
  printf ("          dx = % 15.8f % 15.8f % 15.8f \n",
          gsl_vector_get (s->dx, 0),
          gsl_vector_get (s->dx, 1),
          gsl_vector_get (s->dx, 2));
}
// for two parameters:
void print_state_2 (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3u x = % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          (unsigned int)iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_blas_dnrm2 (s->f));
  printf ("          dx = % 15.8f % 15.8f \n",
          gsl_vector_get (s->dx, 0),
          gsl_vector_get (s->dx, 1));
}

// Fit a Gaussian to provided data
bool StopPow::fit_Gaussian(std::vector<double> & data_x, 
                            std::vector<double> & data_y, 
                            std::vector<double> & data_std, 
                            std::vector<double> & fit,
                            std::vector<double> & fit_unc,
                            double & chi2_dof,
                            bool verbose)
{
	bool ret = true;
    int status;
    unsigned int iter = 0;
    const size_t n = data_x.size(); // number of data points
    const size_t p = 3; // number of parameters

    // Find max point:
    int max_i = find_max_i(data_y);
    double scale = data_y[max_i];
    // GSL becomes unhappy if values being fit are very large.
    // Therefore create scaled data (values of order unity):
    std::vector<double> data_y2(data_y);
    std::vector<double> data_std2(data_std);
    for(int i=0; i<data_y.size(); i++)
    {
        data_y2[i] *= 1./scale;
        data_std2[i] *= 1./scale;
    }

    // some automated guesses to start from:
    double x_init[3] = { data_y[max_i]/scale, data_x[max_i], data_x[1]-data_x[0] };
    // allocate memory for x (fit) values:
    gsl_vector_view x = gsl_vector_view_array (x_init, p);

    // allocate covariance matrix, and set up data struct:
    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    struct data d = { n, data_x, data_y2, data_std2};

    // Set up function for GSL:
    gsl_multifit_function_fdf f;
    f.f = &gauss_f;
    f.df = &gauss_df;
    f.fdf = &gauss_fdf;
    f.n = n;
    f.p = p;
    f.params = &d;

    // Set up the solver:
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc (T, n, p);
    gsl_multifit_fdfsolver_set (s, &f, &x.vector);

    // Iterative loop for the fit:
    do
    {
    	// advance the solution
        iter++;
        status = gsl_multifit_fdfsolver_iterate (s);

        if(verbose)
        {
            printf ("status = %s\n", gsl_strerror (status));
            print_state (iter, s);
        }

        // detect an error:
        if (status)
        {
        	ret = false;
        	break;
        }

        // check if we need to continue
        status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
    }
    while (status == GSL_CONTINUE && iter < 1000);

    if(iter >= 1000)
    	ret = false; // did not converge!

    // Get the covariance matrix and calculate chi^2:
    gsl_multifit_covar (s->J, 1e-4, covar);
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;

    // output if requested:
    if(verbose)
    {
        printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
        printf ("A      = %.5f +/- %.5f\n", gsl_vector_get(s->x,0)*scale, sqrt(gsl_matrix_get(covar,0,0))*scale );
        printf ("mu     = %.5f +/- %.5f\n", gsl_vector_get(s->x,1), sqrt(gsl_matrix_get(covar,1,1)) );
        printf ("sigma  = %.5f +/- %.5f\n", gsl_vector_get(s->x,2), sqrt(gsl_matrix_get(covar,2,2)) );
        printf ("status = %s\n", gsl_strerror (status));
    }

    // Provide results via referenced variables:
    fit.resize(3);
    fit_unc.resize(3);
    chi2_dof = pow(chi, 2.0) / dof;
    for(int j=0; j<3; j++)
    {
    	fit[j] = gsl_vector_get(s->x, j);
    	fit_unc[j] = sqrt(gsl_matrix_get(covar,j,j));
    }
    // correct for scaling of amplitude:
    fit[0] *= scale;
    fit_unc[0] *= scale;

    // free memory:
    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);

    return ret;
}

// Basic fitting to infer rhoR
bool StopPow::fit_rhoR(std::vector<double> & data_x, 
                        std::vector<double> & data_y, 
                        std::vector<double> & data_std,
                        double & dE, 
                        std::vector<double> & fit,
                        std::vector<double> & fit_unc,
                        double & chi2_dof,
                        StopPow & s,
                        double & E0,
                        double & rhoR,
                        double & rhoR_unc,
                        bool verbose)
{
	bool ret = true;

	// Wrap calculation in try/catch
	try
	{
		// Call the Gaussian fitting function:
		ret = fit_Gaussian(data_x, data_y, data_std, fit, fit_unc, chi2_dof, verbose);

		// Negative uncertainty in energy is a bit silly:
		dE = fmax(dE, 0.);

		// Mean energy and total uncertainty
		double E = fit[1];
		double E_unc = sqrt( pow(fit_unc[1],2) + pow(dE,2) ); // add in quadrature

		// get current mode of the model
		int currMode = s.get_mode();
		s.set_mode(s.MODE_RHOR);

		rhoR = s.Thickness(E0, E);
		// calculate and store rhoR error bar:
		double rhoR_min = s.Thickness(E0, E+E_unc);
		double rhoR_max = s.Thickness(E0, E-E_unc);
		rhoR_unc = (rhoR_max - rhoR_min) / 2.;

		if(verbose)
		{
			std::cout << "Fit E = " << E << " +/- " << E_unc << std::endl;
			std::cout << "rhoR = " << rhoR << " +/- " << rhoR_unc << std::endl;
		}

		// reset mode of s:
		s.set_mode(currMode);
	}
	catch(...)
	{
		// something bad happened
		ret = false;
	}

	return ret;
}

// ----------------------------------------------------
//				Forward-Fitting Routine
// ----------------------------------------------------
// data struct for forward fitting
struct ff_data {
  size_t n;
  std::vector<double> & x;
  std::vector<double> & y;
  std::vector<double> & sigma;
  double E0;
  StopPow::StopPow * s;
};

// Gaussian forward-fit function, in form to be used with gsl multifit library
int forward_gauss_f(const gsl_vector * p, void * data, gsl_vector * f) 
{
    // retrieve data from struct passed to this function:
    size_t n = ((struct ff_data *)data)->n;
    std::vector<double> x = ((struct ff_data *)data)->x;
    std::vector<double> y = ((struct ff_data *)data)->y;
    std::vector<double> stdev = ((struct ff_data *) data)->sigma;
    double E0 = ((struct ff_data *) data)->E0;
    StopPow::StopPow * s = ((struct ff_data *) data)->s;

    // fitting parameters:
    double rhoR = gsl_vector_get (p, 0);
    double A = gsl_vector_get(p, 1);
    double sigma = gsl_vector_get(p, 2);

    // loop over data, putting (y_ff[i] - y[i])/sigma[i] into f (i.e. chi for each point)
    for(size_t i=0; i < n; i++)
    {
    	double Ein = s->Ein(x[i], rhoR);
    	double y_eval = (A/(sqrt(2*M_PI)*sigma)) * exp(-1.*pow(Ein-E0,2)/(2*pow(sigma,2)));
    	// correct  for "accordion" effect on spectrum:
    	double accordion = 0.1 / (s->Ein(x[i]+0.05,rhoR) - s->Ein(x[i]-0.05,rhoR));
    	y_eval *= 1./accordion;
    	// store result:
    	gsl_vector_set(f, i, (y_eval-y[i])/stdev[i]);
    }

    return GSL_SUCCESS;
}

// stuff needed to calculate derivatives in following function
struct ff_deriv_params {
	gsl_vector * p2;
	gsl_vector * f;
	std::vector<double> & x;
	std::vector<double> & y;
	std::vector<double> & stdev;
	double E0;
	StopPow::StopPow * s;
	int i;
};

// Jacobian matrix for previous
int forward_gauss_df(const gsl_vector * p, void * data, gsl_matrix * J) 
{
    // retrieve data from params passed to this function
    size_t n = ((struct ff_data *)data)->n;
    std::vector<double> x = ((struct ff_data *)data)->x;
    std::vector<double> y = ((struct ff_data *)data)->y;
    std::vector<double> stdev = ((struct ff_data *) data)->sigma;
    double E0 = ((struct ff_data *) data)->E0;
    StopPow::StopPow * s = ((struct ff_data *) data)->s;

    // fitting parameters:
    double rhoR = gsl_vector_get (p, 0);
    double A = gsl_vector_get(p, 1);
    double sigma = gsl_vector_get(p, 2);

    // allocate some memory (used for calculating derivatives)
	gsl_vector * p2 = gsl_vector_alloc(3);
	gsl_vector_memcpy(p2, p);
	gsl_vector * f = gsl_vector_alloc(1);

	// for calculating derivatives, need these parameters:
	ff_deriv_params params = {p2, f, x, y, stdev, E0, s, 0};

	// Function to be used to calculate derivative with respect to rhoR parameter
	// This whole approach is clunky but robust
	auto drhoR_func = [] (double rhoR, void * params) {
		// cast the parameter struct:
		ff_deriv_params * p = ((ff_deriv_params*) params);
		// rhoR is the free parameter, set it in p2 which will be fed to forward_gauss_f
		gsl_vector_set(p->p2, 0, rhoR); 
		// size one data array
		std::vector<double> x {p->x[p->i]};
		std::vector<double> y {p->y[p->i]};
		std::vector<double> stdev {p->stdev[p->i]};
		ff_data data2 = {1, x, y, stdev, p->E0, p->s};
		forward_gauss_f(p->p2, &data2, p->f);
		return (double)gsl_vector_get(p->f, 0);
	};
	gsl_function deriv_F_rhoR;
	deriv_F_rhoR.function = drhoR_func;
	deriv_F_rhoR.params = &params;

	// Function to be used to calculate derivative with respect to A parameter
	auto dA_func = [] (double A, void * params) {
		// cast the parameter struct:
		ff_deriv_params * p = ((ff_deriv_params*) params);
		// x is the free parameter, set it in p2 which will be fed to forward_gauss_f
		gsl_vector_set(p->p2, 1, A); 
		// size one data array
		std::vector<double> x {p->x[p->i]};
		std::vector<double> y {p->y[p->i]};
		std::vector<double> stdev {p->stdev[p->i]};
		ff_data data2 = {1, x, y, stdev, p->E0, p->s};
		forward_gauss_f(p->p2, &data2, p->f);
		return (double)gsl_vector_get(p->f, 0);
	};
	gsl_function deriv_F_A;
	deriv_F_A.function = dA_func;
	deriv_F_A.params = &params;

	// Function to be used to calculate derivative with respect to sigma parameter
	auto ds_func = [] (double s, void * params) {
		// cast the parameter struct:
		ff_deriv_params * p = ((ff_deriv_params*) params);
		// x is the free parameter, set it in p2 which will be fed to forward_gauss_f
		gsl_vector_set(p->p2, 2, s); 
		// size one data array
		std::vector<double> x {p->x[p->i]};
		std::vector<double> y {p->y[p->i]};
		std::vector<double> stdev {p->stdev[p->i]};
		ff_data data2 = {1, x, y, stdev, p->E0, p->s};
		forward_gauss_f(p->p2, &data2, p->f);
		return (double)gsl_vector_get(p->f, 0);
	};
	gsl_function deriv_F_s;
	deriv_F_s.function = ds_func;
	deriv_F_s.params = &params;

    // calculate derivative for for all points in Jacobian
    for(size_t i=0; i < n; i++)
    {
    	params.i = i;

    	// Reset parameters and calculate df/drhoR
		gsl_vector_memcpy(p2, p);
    	double rhoR_deriv, err;
    	gsl_deriv_central(&deriv_F_rhoR, rhoR, 0.01, &rhoR_deriv, &err);

    	// Reset parameters and calculate df/dA
    	gsl_vector_memcpy(p2, p);
    	double A_deriv;
    	gsl_deriv_central(&deriv_F_A, A, 1e3, &A_deriv, &err);

    	// Reset parameters and calculate df/dsigma
		gsl_vector_memcpy(p2, p);
    	double s_deriv;
    	gsl_deriv_central(&deriv_F_s, sigma, 1e-2, &s_deriv, &err);

        // Set Jacobian:
        gsl_matrix_set (J, i, 0, rhoR_deriv);
        gsl_matrix_set (J, i, 1, A_deriv);
        gsl_matrix_set (J, i, 2, s_deriv);
    }

    gsl_vector_free(p2);
    gsl_vector_free(f);
    return GSL_SUCCESS;
}

// Convenient combination function for above, which GSL needs
int forward_gauss_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
  forward_gauss_f (x, data, f);
  forward_gauss_df (x, data, J);

  return GSL_SUCCESS;
}

// Forward fit a Gaussian to data to infer rhoR
bool StopPow::forward_fit_rhoR(std::vector<double> & data_x, 
                                std::vector<double> & data_y, 
                                std::vector<double> & data_std,
                                double & dE, 
                                double & chi2_dof,
                                StopPow & s,
                                double & E0,
                                std::vector<double> & fit,
                                std::vector<double> & fit_unc,
                                bool verbose)
{
	bool ret = true;
    int status;
    unsigned int iter = 0;
    const size_t n = data_x.size(); // number of data points
    const size_t p = 3; // number of parameters: rhoR, A, sigma
    double chi, dof;
    int mode_init = s.get_mode();
    s.set_mode(s.MODE_RHOR);

    // allocate solver:
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *solver;
    T = gsl_multifit_fdfsolver_lmsder;

    // allocate covariance matrix
    gsl_matrix *covar = gsl_matrix_alloc (p, p);

    // Find max point:
    int max_i = find_max_i(data_y);
    double scale = data_y[max_i];
    // GSL becomes unhappy if values being fit are very large.
    // Therefore create scaled data (values of order unity):
    std::vector<double> data_y2(data_y);
    std::vector<double> data_std2(data_std);
    for(int i=0; i<data_y.size(); i++)
    {
        data_y2[i] *= 1./scale;
        data_std2[i] *= 1./scale;
    }

    // Use a standard Gaussian fit rhoR as an initial guess
    std::vector<double> dummy_fit, dummy_fit_unc;
    double dummy_chi2_dof, dummy_rhoR, dummy_rhoR_unc;
    fit_rhoR(data_x, data_y2, data_std2, dE, dummy_fit, dummy_fit_unc, dummy_chi2_dof, s, E0, dummy_rhoR, dummy_rhoR_unc, false);
    // initial guess:
    double x_init[3] = { dummy_rhoR, dummy_fit[0], dummy_fit[2] };

    // Run the fit routine three times for initial energy, including provided error bar:
    std::vector<double> results;
    std::vector<double> Eshift {-dE, +dE, 0};
    for(auto & fit_dE : Eshift)
    {
    	// for error analysis, shift energies:
    	std::vector<double> data_x2(data_x);
    	for(int i=0; i<n; i++)
    		data_x2[i] += fit_dE;

	    // set up the data for fitting:
	    struct ff_data d = {n, data_x2, data_y2, data_std2, E0, &s};

	    // Set up function for GSL:
	    gsl_multifit_function_fdf f;
	    f.f = &forward_gauss_f;
	    f.df = &forward_gauss_df;
	    f.fdf = &forward_gauss_fdf;
	    f.n = n;
	    f.p = p;
	    f.params = &d;

	    // allocate memory for x (fit) values:
	    gsl_vector_view x = gsl_vector_view_array (x_init, p);

	    // Set up the solver:
	    solver = gsl_multifit_fdfsolver_alloc (T, n, p);
	    gsl_multifit_fdfsolver_set (solver, &f, &x.vector);

	    // Iterative loop for the fit:
        iter = 0;
	    do
	    {
	        iter++;
	        status = gsl_multifit_fdfsolver_iterate (solver);

	        if(verbose)
	        {
	            printf ("status = %s\n", gsl_strerror (status));
	            print_state (iter, solver);
	        }

	        // detect an error:
	        if (status)
	        {
	        	ret = false;
	        	//break;
	        }

	        status = gsl_multifit_test_delta (solver->dx, solver->x, 1e-4, 1e-4);
	    }
	    while (status == GSL_CONTINUE && iter < 100);

	    if(iter >= 100)
	    	ret = false; // did not converge!

	    // Get the covariance matrix and calculate chi^2:
	    gsl_multifit_covar (solver->J, 1e-4, covar);
	    chi = gsl_blas_dnrm2(solver->f);
	    dof = n - p;

	    results.push_back( gsl_vector_get(solver->x,0) );
	}

    // output if requested:
    if(verbose)
    {
        printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
        printf ("rhoR   = %.5f +/- %.5f\n", gsl_vector_get(solver->x,0), sqrt(gsl_matrix_get(covar,0,0)) );
        printf ("A      = %.5f +/- %.5f\n", gsl_vector_get(solver->x,1)*scale, sqrt(gsl_matrix_get(covar,1,1))*scale );
        printf ("sigma  = %.5f +/- %.5f\n", gsl_vector_get(solver->x,2), sqrt(gsl_matrix_get(covar,2,2)) );
        printf ("status = %s\n", gsl_strerror (status));
    }

    // set results:
    fit.reserve(3);
    fit_unc.reserve(3);
    for(int i=0; i<3; i++)
    {
        fit[i] = gsl_vector_get(solver->x, i);
        fit_unc[i] = sqrt(gsl_matrix_get(covar,i,i));
    }
    // rhoR has some extra uncertainty due to dE uncertainty in addition to intrinsic fit unc:
    double drhoR_1 = fabs(results[1]-results[0])/2.;
    double rhoR_unc = sqrt( pow(drhoR_1,2) + gsl_matrix_get(covar,0,0) );
    fit_unc[0] = rhoR_unc;
    // Fix scale for amplitude:
    fit[1] *= scale;
    fit_unc[1] *= scale;

    // free memory:
    gsl_multifit_fdfsolver_free (solver);
    gsl_matrix_free (covar);

    // return s to original state:
    s.set_mode(mode_init);

    return ret;
}


// ----------------------------------------------------
//				Deconvolution Routine
// ----------------------------------------------------
// data structure for deconvolution algorithms
struct dc_data {
  std::vector<double> & x;
  std::vector<double> & y;
  std::vector<double> & sigma;
  double E0;
  StopPow::StopPow * s;
  double fit_unc;
}; 

// Function to be minimized
double deconvolve_f(double rhoR, void * params)
{
	// cast the parameter struct:
	dc_data * p = ((dc_data*) params);

	// Calculate deconvolved spectrum:
	std::vector<double> x(p->x);
	std::vector<double> y(p->y);
	std::vector<double> sigma(p->sigma);
	StopPow::shift(*(p->s), -rhoR, x, y, sigma);

	// Gaussian fit the deconvolved spectrum:
	std::vector<double> fit;
	std::vector<double> fit_unc;
	double chi2;
	StopPow::fit_Gaussian(x, y, sigma, fit, fit_unc, chi2, false);

	// store the uncertainty from Gaussian fit:
	p->fit_unc = fit_unc[1];

	return fit[1] - p->E0;
}

// Use deconvolution to fit/analyze rhoR
bool StopPow::deconvolve_fit_rhoR(std::vector<double> & data_x, 
                                    std::vector<double> & data_y, 
                                    std::vector<double> & data_std,
                                    double & dE, 
                                    double & chi2_dof,
                                    StopPow & s,
                                    double & E0,
                                    std::vector<double> & fit,
                                    std::vector<double> & fit_unc,
                                    bool verbose)
{
	// make sure mode is set to rhoR
    int mode_init = s.get_mode();
    s.set_mode(s.MODE_RHOR);

	int status;
	int iter = 0, max_iter = 100;
	// Set up stuff for GSL root finding:
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *solver;
	gsl_function F;
	F.function = &deconvolve_f;
	T = gsl_root_fsolver_brent;
	solver = gsl_root_fsolver_alloc (T);

	// for iteration:
	double r, x_lo, x_hi;
    // Run the fit routine three times for initial energy, including provided error bar:
    std::vector<double> results;
    std::vector<double> results_unc;
    std::vector<double> Eshift {-dE, +dE, 0};
    for(auto & fit_dE : Eshift)
    {
    	// for error analysis, shift energies:
    	std::vector<double> data_x2(data_x);
    	for(int i=0; i<data_x2.size(); i++)
    		data_x2[i] += fit_dE;

    	// to feed into fitting function:
		struct dc_data params = {data_x2, data_y, data_std, E0, &s, 0};
		F.params = &params;

		// get an initial guess:
		std::vector<double> fit, fit_unc;
		double guess, guess_unc;
		fit_rhoR(data_x2, data_y, data_std, dE, fit, fit_unc, chi2_dof, s, E0, guess, guess_unc, false);
		// set the solver initial point
		gsl_root_fsolver_set (solver, &F, guess*0.75, 1.25*guess);

		if(verbose)
		{
			printf ("using %s method\n", 
			      gsl_root_fsolver_name (solver));

			printf ("%5s [%9s, %9s] %9s %10s %9s\n",
			      "iter", "lower", "upper", "root", 
			      "err", "err(est)");
		}

        iter = 0;
		do
		{
			// iterate the solver:
			iter++;
			status = gsl_root_fsolver_iterate (solver);
			r = gsl_root_fsolver_root (solver);
			x_lo = gsl_root_fsolver_x_lower (solver);
			x_hi = gsl_root_fsolver_x_upper (solver);

			// convergence check
			status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);

			if(verbose)
			{
				if (status == GSL_SUCCESS)
					printf ("Converged:\n");

				printf ("%5d [%.7f, %.7f] %+.7f %.7f\n",
				      iter, x_lo, x_hi,
				      r, 
				      x_hi - x_lo);
			}
		}
		while (status == GSL_CONTINUE && iter < max_iter);

		// store results
		results.push_back(gsl_root_fsolver_root(solver));
		results_unc.push_back(params.fit_unc);
	}

	// set final results:
	double rhoR = results[2];
    // error in rhoR is combination of uncertainty due to dE, intrinsic fit unc, and root-finding convergence
    double drhoR_1 = fabs(results[1]-results[0])/2.;
    double rhoR_unc = sqrt( pow(drhoR_1,2) + pow(x_hi-x_lo,2) + pow(results_unc[2],2) );

    // Also run a Gaussian fit to get (explicitly) the chi2 and other fit parameters:
    std::vector<double> data_x2(data_x);
    std::vector<double> data_y2(data_y);
    std::vector<double> data_sigma2(data_std);
    shift(s, -rhoR, data_x2, data_y2, data_sigma2);

    // Gaussian fit the deconvolved spectrum:
    fit_Gaussian(data_x2, data_y2, data_sigma2, fit, fit_unc, chi2_dof, false);

    // Have to re-jigger fit. Gaussian above puts [A,mu,sigma], we want [rhoR,A,sigma]:
    fit[1] = fit[0];
    fit_unc[1] = fit_unc[0];
    fit[0] = rhoR;
    fit_unc[0] = rhoR_unc;

    // free memory
	gsl_root_fsolver_free (solver);

	// restore mode of dE/dx model:
	s.set_mode(mode_init);

	return (status == GSL_SUCCESS);
}


// ----------------------------------------------------
//              Stopping power fit
// ----------------------------------------------------
// data struct for forward fitting to dE/dx
struct ff_dEdx_data {
  std::vector<double> & x;
  std::vector<double> & y;
  std::vector<double> & sigma;
  double E0;
  double sigma0;
  double rhoR;
  StopPow::StopPow_Fit * s;
};

// Gaussian forward-fit function, in form to be used with gsl multifit library
int forward_dEdx_f(const gsl_vector * p, void * data, gsl_vector * f) 
{
    // retrieve data from struct passed to this function:
    ff_dEdx_data * d = ((struct ff_dEdx_data *)data);
    // std::vector<double> x = d->x;
    // std::vector<double> y = d->y;
    // std::vector<double> stdev = d->sigma;
    // double E0 = d->E0;
    // double sigma0 = d->sigma0;
    // double rhoR = d->rhoR;
    // StopPow::StopPow_Fit * s = d->s;

    // fitting parameters:
    double factor = gsl_vector_get (p, 0);
    double A = gsl_vector_get(p, 1);
    // use factor to set s:
    d->s->set_factor(factor);

    // loop over data, putting (y_ff[i] - y[i])/sigma[i] into f (i.e. chi for each point)
    double Ein, y_eval, accordion;
    for(size_t i=0; i < d->x.size(); i++)
    {
        Ein = d->s->Ein(d->x[i], d->rhoR);
        y_eval = (A/(sqrt(2*M_PI)*d->sigma0)) * exp(-1.*pow(Ein-d->E0,2)/(2*pow(d->sigma0,2)));
        // correct  for "accordion" effect on spectrum:
        accordion = 0.1 / (d->s->Ein(d->x[i]+0.05,d->rhoR) - d->s->Ein(d->x[i]-0.05,d->rhoR));
        y_eval *= 1./accordion;
        // store result:
        gsl_vector_set(f, i, (y_eval - d->y[i]) / d->sigma[i]);
    }

    return GSL_SUCCESS;
}

// stuff needed to calculate derivatives in following function
struct ff_dEdx_deriv_params {
    ff_dEdx_data * d;
    gsl_vector * p2;
    gsl_vector * f;
    int i;
    ff_dEdx_data * data2;
};

// Jacobian matrix for previous
int forward_dEdx_df(const gsl_vector * p, void * data, gsl_matrix * J) 
{
    // retrieve data from struct passed to this function:
    ff_dEdx_data * d = ((struct ff_dEdx_data *)data);

    // fitting parameters:
    double factor = gsl_vector_get (p, 0);
    double A = gsl_vector_get(p, 1);

    // allocate some memory (used for calculating derivatives)
    gsl_vector * p2 = gsl_vector_alloc(2);
    gsl_vector_memcpy(p2, p);
    gsl_vector * f = gsl_vector_alloc(1);
    std::vector<double> x {d->x[0]};
    std::vector<double> y {d->y[0]};
    std::vector<double> s {d->sigma[0]};
    ff_dEdx_data data2 = {x, y, s, d->E0, d->sigma0, d->rhoR, d->s};

    // for calculating derivatives, need these parameters:
    ff_dEdx_deriv_params params = {d, p2, f, 0, &data2};

    // Function to be used to calculate derivative with respect to 'factor' parameter
    // This whole approach is clunky but robust
    auto dfactor_func = [] (double factor, void * params) {
        // cast the parameter struct:
        ff_dEdx_deriv_params * p = ((ff_dEdx_deriv_params*) params);
        // factor is the free parameter, set it in p2 which will be fed to forward_dEdx_f
        gsl_vector_set(p->p2, 0, factor); 
        // size one data array
        // make new data structure to be passed in calculating f at one point (corresponding to index i)
        int i = p->i;
        p->data2->x[0] = p->d->x[i];
        p->data2->y[0] = p->d->y[i];
        p->data2->sigma[0] = p->d->sigma[i];
        forward_dEdx_f(p->p2, p->data2, p->f);
        return (double)gsl_vector_get(p->f, 0);
    };
    gsl_function deriv_F_factor;
    deriv_F_factor.function = dfactor_func;
    deriv_F_factor.params = &params;

    // Function to be used to calculate derivative with respect to A parameter
    auto dA_func = [] (double A, void * params) {
        // cast the parameter struct:
        ff_dEdx_deriv_params * p = ((ff_dEdx_deriv_params*) params);
        // A is the free parameter, set it in p2 which will be fed to forward_dEdx_f
        gsl_vector_set(p->p2, 1, A); 
        // size one data array
        // make new data structure to be passed in calculating f at one point (corresponding to index i)
        int i = p->i;
        p->data2->x[0] = p->d->x[i];
        p->data2->y[0] = p->d->y[i];
        p->data2->sigma[0] = p->d->sigma[i];
        forward_dEdx_f(p->p2, p->data2, p->f);
        return (double)gsl_vector_get(p->f, 0);
    };
    gsl_function deriv_F_A;
    deriv_F_A.function = dA_func;
    deriv_F_A.params = &params;

    // calculate derivative for for all points in Jacobian
    for(size_t i=0; i < d->x.size(); i++)
    {
        params.i = i;

        // Reset parameters and calculate df/drhoR
        gsl_vector_memcpy(p2, p);
        double factor_deriv, err;
        gsl_deriv_central(&deriv_F_factor, factor, 0.01, &factor_deriv, &err);

        // Reset parameters and calculate df/dA
        gsl_vector_memcpy(p2, p);
        double A_deriv;
        gsl_deriv_central(&deriv_F_A, A, 1e3, &A_deriv, &err);

        // Set Jacobian:
        gsl_matrix_set (J, i, 0, factor_deriv);
        gsl_matrix_set (J, i, 1, A_deriv);
    }

    gsl_vector_free(p2);
    gsl_vector_free(f);
    return GSL_SUCCESS;
}

// Convenient combination function for above, which GSL needs
int forward_dEdx_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
  forward_dEdx_f (x, data, f);
  forward_dEdx_df (x, data, J);

  return GSL_SUCCESS;
}

// Main method, which actually does the fitting routine
bool StopPow::forward_fit_dEdx(std::vector<double> & data_x, 
                                std::vector<double> & data_y, 
                                std::vector<double> & data_std,
                                double & dE, 
                                double E0,
                                double E0_unc,
                                double sigma,
                                double sigma_unc,
                                double rhoR,
                                double rhoR_unc,
                                StopPow_Fit & s,
                                double & chi2_dof,
                                std::vector<double> & fit,
                                std::vector<double> & fit_unc,
                                bool verbose)
{

    bool ret = true;
    int status;
    unsigned int iter = 0;
    const size_t n = data_x.size(); // number of data points
    const size_t p = 2; // number of parameters: factor, A
    double chi, dof;
    // set mode to rhoR
    int mode_init = s.get_mode();
    s.set_mode(s.MODE_RHOR);

    // allocate solver:
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *solver;
    T = gsl_multifit_fdfsolver_lmsder;

    // allocate covariance matrix
    gsl_matrix *covar = gsl_matrix_alloc (p, p);

    // Find max point:
    int max_i = find_max_i(data_y);
    double scale = data_y[max_i];
    // GSL becomes unhappy if values being fit are very large.
    // Therefore create scaled data (values of order unity):
    std::vector<double> data_y2(data_y);
    std::vector<double> data_std2(data_std);
    for(int i=0; i<data_y.size(); i++)
    {
        data_y2[i] *= 1./scale;
        data_std2[i] *= 1./scale;
    }

    // Use a standard Gaussian fit as an initial guess for amplitude (yield)
    std::vector<double> dummy_fit, dummy_fit_unc;
    double dummy_chi2;
    fit_Gaussian(data_x, data_y2, data_std2, dummy_fit, dummy_fit_unc, dummy_chi2, false);
    // initial guess. Factor should be ~ 1
    double x_init[2] = { 1, dummy_fit[0] };

    // Run the fit routine five times, which allow varying:
    // E of data (via dE uncertainty)
    // E0 initial energy
    // sigma0 initial width
    // nominal case is last one for convenience in using covar matrix
    std::vector<double> results;
    std::vector<double> vary_E {-dE, +dE, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> vary_E0 {E0, E0, E0-E0_unc, E0+E0_unc, E0, E0, E0, E0, E0};
    std::vector<double> vary_sigma0 {sigma, sigma, sigma, sigma, sigma-sigma_unc, sigma+sigma_unc, sigma, sigma, sigma}; 
    std::vector<double> vary_rhoR {rhoR, rhoR, rhoR, rhoR, rhoR, rhoR, rhoR+rhoR_unc, rhoR-rhoR_unc, rhoR};
    for(int i=0; i<9; i++) // loop corresponds to varying above params
    {
        // for error analysis, shift energies:
        std::vector<double> data_x2(data_x);
        for(int j=0; j<n; j++)
            data_x2[j] += vary_E[i];

        // set up the data for fitting:
        struct ff_dEdx_data d = {data_x2, data_y2, data_std2, vary_E0[i], vary_sigma0[i], vary_rhoR[i], &s};

        // Set up function for GSL:
        gsl_multifit_function_fdf f;
        f.f = &forward_dEdx_f;
        f.df = &forward_dEdx_df;
        f.fdf = &forward_dEdx_fdf;
        f.n = n;
        f.p = p;
        f.params = &d;

        // allocate memory for x (fit) values:
        gsl_vector_view x = gsl_vector_view_array (x_init, p);

        // Set up the solver:
        solver = gsl_multifit_fdfsolver_alloc (T, n, p);
        gsl_multifit_fdfsolver_set (solver, &f, &x.vector);

        // Iterative loop for the fit:
        iter = 0;
        do
        {
            iter++;
            status = gsl_multifit_fdfsolver_iterate (solver);

            if(verbose)
            {
                printf ("status = %s\n", gsl_strerror (status));
                print_state_2 (iter, solver);
            }

            // detect an error:
            if (status)
            {
                ret = false;
                break;
            }

            status = gsl_multifit_test_delta (solver->dx, solver->x, 1e-3, 1e-3);
        }
        while (status == GSL_CONTINUE && iter < 100);

        if(iter >= 100)
            ret = false; // did not converge!

        // Get the covariance matrix and calculate chi^2:
        gsl_multifit_covar (solver->J, 1e-3, covar);
        chi = gsl_blas_dnrm2(solver->f);
        dof = n - p;

        results.push_back( gsl_vector_get(solver->x,0) );
    }

    // output if requested:
    if(verbose)
    {
        printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
        printf ("factor   = %.5f +/- %.5f\n", gsl_vector_get(solver->x,0), sqrt(gsl_matrix_get(covar,0,0)) );
        printf ("A      = %.5f +/- %.5f\n", gsl_vector_get(solver->x,1)*scale, sqrt(gsl_matrix_get(covar,1,1))*scale );
        printf ("status = %s\n", gsl_strerror (status));
    }

    // set results:
    fit.clear();
    fit.push_back( gsl_vector_get(solver->x, 0) );
    fit.push_back( gsl_vector_get(solver->x, 1)*scale );
    // error in rhoR is combination of uncertainty due to dE and intrinsic fit unc:
    fit_unc.clear();
    double dfactor_1 = fabs(results[1]-results[0])/2.;
    double dfactor_2 = fabs(results[3]-results[2])/2.;
    double dfactor_3 = fabs(results[5]-results[4])/2.;
    double dfactor_4 = fabs(results[7]-results[6])/2.;
    fit_unc.push_back( sqrt( pow(dfactor_1,2) + pow(dfactor_2,2) + pow(dfactor_3,2) + pow(dfactor_4,2) + gsl_matrix_get(covar,0,0) ) );
    fit_unc.push_back( sqrt(gsl_matrix_get(covar,1,1))*scale );

    // free memory:
    gsl_multifit_fdfsolver_free (solver);
    gsl_matrix_free (covar);

    // return s to original state:
    s.set_mode(mode_init);

    return ret;
}
