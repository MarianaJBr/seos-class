#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <stdlib.h>

typedef struct
{
    double wa_fld;
    double q_fld;
    double zt_fld;
} background_w_fl_i_args;

int background_w_fld(
    double a0,
    double w0,
    double wa,
    double q,
    double zt,
    double a,
    double a_today,
    double *w_fld,
    double *dw_over_da_fld,
    double *integral_fld)
{

    /** - first, define the function w(a) */
    /** mjb:seos: here we define seos instead of CPL*/
    /**w_fld = pba->w0_fld + pba->wa_fld * (1. - a / pba->a_today);
     */
    /** mjb proper way to introduce local variables*/
    // double a0 = pba->a_today;
    // double w0 = pba->w0_fld;
    // double wa = pba->wa_fld;
    // double q = pba->q_fld;
    // double zt = pba->zt_fld;

    /**w_fld =  w0 + wa * pow((a0 - a), q) / (pow(a * zt, q) + pow((a0 - a), q));*/

    /* mjb: seos: add conditional to avoid divergence */
    if (a > a_today)
    {
        *w_fld = w0;
    }
    else
    {
        *w_fld = w0 + wa * pow((a0 - a), q) / (pow(a * zt, q) + pow((a0 - a), q));
    }

    /** - then, give the corresponding analytic derivative dw/da (used
        by perturbation equations; we could compute it numerically,
        but with a loss of precision; as long as there is a simple
        analytic expression of the derivative of the previous
        function, let's use it! */
    /*mjb:seos: we add the corresponding derivative of the eos w(a) wrt a */
    /* numerator = a0*wa*q*(a*zt)^q*(a0-a)^(q-1))  */
    /* denominator = a ( (a0-a)^q + (a*zt)^q )^2  */
    /* dw/da = -numerator/denominator  */
    /**dw_over_da_fld =
            -(a0 * wa * q * (pow(a * zt, q)) * (pow((a0 - a), (q - 1.)))) /
            (a * pow((pow((a * zt), q) + pow((a0 - a), q)), 2.));*/

    /* mjb:seos: add conditional to avoid divergence. TODO: necessary ?*/
    if (a > a_today)
    {
        *dw_over_da_fld = 0.;
    }
    else
    {
        *dw_over_da_fld =
            -(a0 * wa * q * (pow(a * zt, q)) * (pow((a0 - a), (q - 1.)))) /
            (a * pow((pow((a * zt), q) + pow((a0 - a), q)), 2.));
    }

    /** - finally, give the analytic solution of the following integral:
        \f$ \int_{a}^{a0} da 3(1+w_{fld})/a \f$. This is used in only
        one place, in the initial conditions for the background, and
        with a=a_ini. If your w(a) does not lead to a simple analytic
        solution of this integral, no worry: instead of writing
        something here, the best would then be to leave it equal to
        zero, and then in background_initial_conditions() you should
        implement a numerical calculation of this integral only for
        a=a_ini, using for instance Romberg integration. It should be
        fast, simple, and accurate enough. */

    /*mjb:seos: following the suggestion we make this equal to zero here*/
    /*
  *integral_fld = 3.*((1.+pba->w0_fld+pba->wa_fld)*log(pba->a_today/a) + pba->wa_fld*(a/pba->a_today-1.));
*/
    *integral_fld = 0.;

    /** note: of course you can generalise these formulas to anything,
      defining new parameters pba->w..._fld. Just remember that so
      far, HyRec explicitly assumes that w(a)= w0 + wa (1-a/a0); but
      Recfast does not assume anything */
    /*mjb:seos: have this in mind?*/

    return EXIT_SUCCESS;
}

double background_w_fl_i(
    double lna,
    void *params)
{

    // void pointer params can point to any type of data
    // We make a type casting
    background_w_fl_i_args *args = (background_w_fl_i_args *)params;

    //mjb: new parameters for eos, q and zt in exponential form:
    double wa = args->wa_fld;
    double q = args->q_fld;
    double zt = args->zt_fld;

    return wa * pow(1 - exp(lna), q) / (pow(zt * exp(lna), q) + pow(1 - exp(lna), q));
}

int main(int argc, char const *argv[])
{
    /* code */
    double integral_seos;
    /* Integrator variables */
    double w0, a;
    double wa_fld;
    double q_fld;
    double zt_fld;
    double integral_fld;
    double result;
    double H0, Omega_fld, rho_fld_today;

    double quad_result, error;
    double eps_abs = 1e-12; //mjb: made this values smaller?//
    double eps_rel = 1e-16;
    size_t max_iter = 1000;


    Omega_fld = 1. - 0.3089;
    H0 = 67.74;
    rho_fld_today = Omega_fld * H0 * H0;

    a = 1e-14;
    w0 = -1.,
    wa_fld = 0.5;
    q_fld = 1;
    zt_fld = 10;
    
    // GSL data structures.
    gsl_integration_workspace *ws = gsl_integration_workspace_alloc(max_iter);

    /*mjb: at this point I need to define the background_w_fl_i*/
    // Here we initialize the integrand arguments.
    background_w_fl_i_args i_args;
    i_args.wa_fld = wa_fld;
    i_args.q_fld = q_fld;
    i_args.zt_fld = zt_fld;

    gsl_function integrand_spec;
    integrand_spec.params = &i_args;
    integrand_spec.function = &background_w_fl_i;

    double val = background_w_fl_i(log(0.1), integrand_spec.params);
    printf("Ival: %,8f\n", val);

    /* integral --> result_integral */
    gsl_integration_qag(&integrand_spec, log(1.0), log(a),
                        eps_abs, eps_rel, max_iter, GSL_INTEG_GAUSS61,
                        ws, &quad_result, &error);

    integral_seos = quad_result;
    // TODO check this sign!
    integral_fld = integral_seos; //  + 3 * (1. + w0) * log(a); //mjb:seos
    /* rho_fld at initial time */
    result = rho_fld_today * exp(integral_fld);

    // Deallocate the used memory.
    gsl_integration_workspace_free(ws);

    printf("rho_fld_today: %.8g\n", rho_fld_today);
    printf("Result: %.8g\n", result);
    printf("Result: %.8g\n", result / rho_fld_today);

    return 0;
}
