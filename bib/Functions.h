#ifndef Functions
#define Functions

double crystalball_function(double x, double mean, double sigma, double alpha, double n)
{
  // evaluate the crystal ball function
  if (sigma < 0.)
    return 0.;
  double z = (x - mean) / sigma;
  if (alpha < 0)
    z = -z;
  double abs_alpha = std::abs(alpha);
  double C = n / abs_alpha * 1. / (n - 1.) * std::exp(-alpha * alpha / 2.);
  double D = std::sqrt(M_PI / 2.) * (1. + ROOT::Math::erf(abs_alpha / std::sqrt(2.)));
  double N = 1. / (sigma * (C + D));
  if (z > -abs_alpha)
    return N * std::exp(-0.5 * z * z);
  else
  {
    // double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
    double nDivAlpha = n / abs_alpha;
    double AA = std::exp(-0.5 * abs_alpha * abs_alpha);
    double B = nDivAlpha - abs_alpha;
    double arg = nDivAlpha / (B - z);
    return N * AA * std::pow(arg, n);
  }
}

double gaussexp_function(double x, double mean, double sigma, double k)
{
  // evaluate the crystal ball function
  if (sigma < 0.)
    return 0.;
  double z = (x - mean) / sigma;
  if (k < 0)
    z = -z;
  double abs_k = std::abs(k);

  if (z >= -abs_k)
    return std::exp(-0.5 * z * z);
  else
  {
    double k2 = 0.5 * abs_k * abs_k;
    return std::exp(k2 + k * z);
  }
}

// Function to compute the decay distribution
double dNdk0(double k0)
{
  double x = k0 / (3.096916 / 2);
  if (k0 < 0.00001)
    return 0;
  if (x > 0.9999)
    return 0;
  return (2 * (1.0 / 137.036) / M_PI) * (2 * log(2 * (3.096916 / 2) / 0.000511 * sqrt(1 - x)) - 1) * (1 / k0) * (1 - x + x * x / 2);
}

// Gaussian function used in convolution
double gaussian(double z, double sigma)
{
  return exp(-z * z / (2 * sigma * sigma)) / (sqrt(2 * M_PI) * sigma);
}

// f(x) as defined in Python code
double f(double x, double br, double eps, double Norm)
{
  if (x < 0.000001)
    return 0.0;
  else if (x >= eps)
    return br * Norm * dNdk0(x);
  else
    return (1 - br) / eps;
}

// F(x) - Convolution of f(x) with Gaussian
double F(double x, double sigma, double br, double eps, double Norm)
{
  double result = 0.0;
  double xmin = eps;
  double xmax = fmax(eps, x + 3.5 * sigma); // Extend to capture the significant part of Gaussian
  if (xmin >= xmax)
    return 0;
  if (x >= 3.0 / 2)
    return 0;

  for (double t = xmin; t <= xmax; t += (xmax - xmin) / 1000)
  {
    result += f(t, br, eps, Norm) * gaussian(x - t, sigma);
  }
  result *= (xmax - xmin) / 1000;

  result += gaussian(x, sigma) * (1 - br);

  // printf("F(x) xmin= %lf\n", xmin);
  // printf("F(x) xmax= %lf\n", xmax);
  // printf("F(x)    x= %lf\n", x);

  return result;
}

// Function to compute convolution for given mee, M, sigma, and branching ratio
double dNdMee(double mee, double M, double sigma, double br)
{
  double eps = 0.001;
  double x = M - mee;
  double Norm = 2.070655621874216; // Normalization constant from the paper

  return F(x, sigma, br, eps, Norm);
}


#endif /* DNDMEE_H */
