# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR

cdef extern from "gsl/gsl_statistics_double.h":
  double gsl_stats_mean ( double data[],  size_t stride,  size_t n)
  double gsl_stats_variance ( double data[],  size_t stride,  size_t n)
  double gsl_stats_sd ( double data[],  size_t stride,  size_t n)
  double gsl_stats_variance_with_fixed_mean ( double data[],  size_t stride,  size_t n,  double mean)
  double gsl_stats_sd_with_fixed_mean ( double data[],  size_t stride,  size_t n,  double mean)
  double gsl_stats_absdev ( double data[],  size_t stride,  size_t n)
  double gsl_stats_skew ( double data[],  size_t stride,  size_t n)
  double gsl_stats_kurtosis ( double data[],  size_t stride,  size_t n)
  double gsl_stats_lag1_autocorrelation ( double data[],  size_t stride,  size_t n)

  double gsl_stats_covariance ( double data1[],  size_t stride1, double data2[],  size_t stride2,  size_t n)

  double gsl_stats_variance_m ( double data[],  size_t stride,  size_t n,  double mean)
  double gsl_stats_sd_m ( double data[],  size_t stride,  size_t n,  double mean)
  double gsl_stats_absdev_m ( double data[],  size_t stride,  size_t n,  double mean)
  double gsl_stats_skew_m_sd ( double data[],  size_t stride,  size_t n,  double mean,  double sd)
  double gsl_stats_kurtosis_m_sd ( double data[],  size_t stride,  size_t n,  double mean,  double sd)
  double gsl_stats_lag1_autocorrelation_m ( double data[],  size_t stride,  size_t n,  double mean)

  double gsl_stats_covariance_m ( double data1[],  size_t stride1, double data2[],  size_t stride2,  size_t n,  double mean1,  double mean2)

  # DEFINED FOR FLOATING POINT TYPES ONLY

  double gsl_stats_wmean ( double w[],  size_t wstride,  double data[],  size_t stride,  size_t n)
  double gsl_stats_wvariance ( double w[],  size_t wstride,  double data[],  size_t stride,  size_t n)
  double gsl_stats_wsd ( double w[],  size_t wstride,  double data[],  size_t stride,  size_t n)
  double gsl_stats_wvariance_with_fixed_mean ( double w[],  size_t wstride,  double data[],  size_t stride,  size_t n,  double mean)
  double gsl_stats_wsd_with_fixed_mean ( double w[],  size_t wstride,  double data[],  size_t stride,  size_t n,  double mean)
  double gsl_stats_wabsdev ( double w[],  size_t wstride,  double data[],  size_t stride,  size_t n)
  double gsl_stats_wskew ( double w[],  size_t wstride,  double data[],  size_t stride,  size_t n)
  double gsl_stats_wkurtosis ( double w[],  size_t wstride,  double data[],  size_t stride,  size_t n)

  double gsl_stats_wvariance_m ( double w[],  size_t wstride,  double data[],  size_t stride,  size_t n,  double wmean)
  double gsl_stats_wsd_m ( double w[],  size_t wstride,  double data[],  size_t stride,  size_t n,  double wmean)
  double gsl_stats_wabsdev_m ( double w[],  size_t wstride,  double data[],  size_t stride,  size_t n,  double wmean)
  double gsl_stats_wskew_m_sd ( double w[],  size_t wstride,  double data[],  size_t stride,  size_t n,  double wmean,  double wsd)
  double gsl_stats_wkurtosis_m_sd ( double w[],  size_t wstride,  double data[],  size_t stride,  size_t n,  double wmean,  double wsd)


  double gsl_stats_pvariance ( double data1[],  size_t stride1,  size_t n1,  double data2[],  size_t stride2,  size_t n2)
  double gsl_stats_ttest ( double data1[],  size_t stride1,  size_t n1,  double data2[],  size_t stride2,  size_t n2)

  double gsl_stats_max ( double data[],  size_t stride,  size_t n)
  double gsl_stats_min ( double data[],  size_t stride,  size_t n)
  void gsl_stats_minmax (double * min, double * max,  double data[],  size_t stride,  size_t n)

  size_t gsl_stats_max_index ( double data[],  size_t stride,  size_t n)
  size_t gsl_stats_min_index ( double data[],  size_t stride,  size_t n)
  void gsl_stats_minmax_index (size_t * min_index, size_t * max_index,  double data[],  size_t stride,  size_t n)

  double gsl_stats_median_from_sorted_data ( double sorted_data[],  size_t stride,  size_t n)
  double gsl_stats_quantile_from_sorted_data ( double sorted_data[],  size_t stride,  size_t n,  double f)

