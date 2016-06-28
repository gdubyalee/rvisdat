library(compiler)
library(inline)

NeutralDrift <- function(lambda, Ns, time.intervals)
{
  theory.timeseries <- matrix(0, Ns+1, length(time.intervals))
  m       <- 1:(Ns-1)
  for(i in 1:length(time.intervals))
  {
    t <- time.intervals[i]
    # For n = 0
    theory.timeseries[1,i] <- 2/Ns*sum((cos(0.5*pi*m/Ns))^2*(1-exp(-4*lambda*t*(sin(0.5*pi*m/Ns))^2)))
    # For n = 1:(Ns-1)
    for(n in 1:(Ns-1))
    {
      theory.timeseries[n+1,i] <- 2/Ns*sum(sin(pi*m/Ns)*sin(pi*m*n/Ns)*exp(-4*lambda*t*(sin(0.5*pi*m/Ns))^2))
    }
    # For n = Ns
    theory.timeseries[Ns+1,i] <- 2/Ns*sum((-1)^(m+1)*(cos(0.5*pi*m/Ns))^2*(1-exp(-4*lambda*t*(sin(0.5*pi*m/Ns))^2)))
  }  
  # Return absolute value (solution to numericaly unstable cases)
  abs(theory.timeseries)
}

CA30Function_MonoClonal <- function(alpha, lambda, Ns, time.intervals)
{
  theory.timeseries <- rep(0, length(time.intervals))
  m       <- 1:(Ns-1)
  intercept <- 0.5*alpha*sum((-1)^(m+1)*(tan(0.5*pi*m/Ns))^(-2))
  aplha_lamb <- alpha*lambda
  for(i in 1:length(time.intervals))
  {
   t <- time.intervals[i]
   theory.timeseries[i] <- 0.5*alpha*sum((-1)^(m+1)*(tan(0.5*pi*m/Ns))^(-2)*(1-exp(-4*lambda*t*(sin(0.5*pi*m/Ns))^2)))
  }  
  aplha_lamb*time.intervals - theory.timeseries
}

CA30Function_MonoClonal_approx <- function(alpha, lambda, Ns, time.intervals)
{
  m          <- 1:(Ns-1)
#   sum((-1)^(m+1)*(tan(0.5*pi*m/Ns))^(-2))
  intercept  <- 0.5*alpha*sum((-1)^(m+1)*(tan(0.5*pi*m/Ns))^(-2))
  aplha_lamb <- alpha*lambda
  # Return value
  print(paste("Slope", aplha_lamb))
  print(paste("Intercept", intercept))
  aplha_lamb*time.intervals - intercept
}


CA30Function_PartialClones <- function(alpha, lambda, Ns, time.intervals)
{
  theory.timeseries   <- rep(0, length(time.intervals))
  theory.timeseries.n <- rep(0, length(time.intervals))
  m       <- 1:(Ns-1)
  n.all   <- 1:(Ns-1)
  for(n in n.all)
  {
    for(i in 1:length(time.intervals))
    {
    t                      <- time.intervals[i]
    theory.timeseries.n[i] <- 0.5*alpha*sum(sin(pi*m/Ns)*sin(pi*m*n/Ns)/(sin(0.5*pi*m/Ns))^2*exp(-4*lambda*t*(sin(0.5*pi*m/Ns))^2))
    }  
  theory.timeseries <- theory.timeseries+theory.timeseries.n   
  }
  # Return value
  alpha*Ns*(Ns-1)/2 - theory.timeseries
}


CA30Function_Intercept_all <- function(alpha, lambda, Ns)
{
  m          <- 1:(Ns-1)
  intercept  <- 0.5*alpha*sum((-1)^(m+1)*(tan(0.5*pi*m/Ns))^(-2))
  alpha*Ns*(Ns-1)/2 - intercept
}



# CA30Adenomas_PartialClones <- function(a, alpha, lambda, Ns, time.intervals)
# {
#   m       <- 1:(Ns-1)
#   n.all   <- 1:(Ns-1)
#   k_0 <- alpha*lambda*Ns
#   cte <- 2*k_0*a/Ns
#   theory.timeseries   <- matrix(0, Ns-1, length(time.intervals))
#   for(n in n.all)
#   {
#     f_nm <- sin(pi*m/Ns)*sin(pi*m*n/Ns)
#     g_nm <- 4*lambda*(sin(0.5*pi*m/Ns))^2
#     for(i in 1:length(time.intervals))
#     {
#       t                       <- time.intervals[i]
#       theory.timeseries[n, i] <- sum(cte*f_nm*(g_nm^(-2) - (t/g_nm + g_nm^(-2))*exp(-g_nm*t)))
#     }  
#   }
#   # Return absolute value (solution to numericaly unstable cases)
#   abs(theory.timeseries)
# }
# 
# CA30Adenomas_FullClones <- function(a, alpha, lambda, Ns, time.intervals)
# {
#   m       <- 1:(Ns-1)
#   k_0     <- alpha*lambda*Ns
#   cte     <- 2*k_0*a/Ns
#   theory.timeseries <- matrix(0, 1, length(time.intervals))
#   f_nm <- (-1)^(m+1)*(cos(0.5*pi*m/Ns))^2
#   g_nm <- 4*lambda*(sin(0.5*pi*m/Ns))^2
#   for(i in 1:length(time.intervals))
#   {
#     t                    <- time.intervals[i]
#     aux1 <- exp(log(t/g_nm + g_nm^(-2))-g_nm*t) 
#     theory.timeseries[i] <- sum(cte*f_nm*(t^2/2 - g_nm^(-2) + aux1))
# #     theory.timeseries[i] <- sum(cte*f_nm*(t^2/2 - g_nm^(-2) + (t/g_nm + g_nm^(-2))*exp(-g_nm*t)))
#   }  
#   # Return absolute value (solution to numericaly unstable cases)
#   abs(theory.timeseries)
# }
# 
# 
# CA30Adenomas_AllClones <- function(a, alpha, lambda, Ns, time.intervals)
# {
#   m       <- 1:(Ns-1)
#   n.all   <- 1:(Ns-1)
#   k_0     <- alpha*lambda*Ns
#   cte     <- 2*k_0*a/Ns
#   theory.timeseries1 <- matrix(0, Ns-1, length(time.intervals))
#   for(n in n.all)
#   {
#     f_nm <- sin(pi*m/Ns)*sin(pi*m*n/Ns)
#     g_nm <- 4*lambda*(sin(0.5*pi*m/Ns))^2
#     for(i in 1:length(time.intervals))
#     {
#       t                       <- time.intervals[i]
#       theory.timeseries1[n, i] <- sum(cte*f_nm*(g_nm^(-2) - (t/g_nm + g_nm^(-2))*exp(-g_nm*t)))
#     }  
#   }
#   # Return value
# #   theory.timeseries1
# 
#   theory.timeseries2 <- matrix(0, 1, length(time.intervals))
#   f_nm <- (-1)^(m+1)*(cos(0.5*pi*m/Ns))^2
#   g_nm <- 4*lambda*(sin(0.5*pi*m/Ns))^2
#   for(i in 1:length(time.intervals))
#   {
#     t                    <- time.intervals[i]
#     theory.timeseries2[i] <- sum(cte*f_nm*(t^2/2 - g_nm^(-2) + (t/g_nm + g_nm^(-2))*exp(-g_nm*t)))
#   }  
#   # Return absolute value (solution to numericaly unstable cases)
#   abs(as.vector(colSums(theory.timeseries1)+theory.timeseries2))
# }
# 
# 
# 
# # CA30Adenomas_AllClones.old.txt <- 
# # "
# #   using namespace arma;
# #   int                Ns_c = as<int>(Ns);
# #   double          alpha_c = as<double>(alpha);
# #   double         lambda_c = as<double>(lambda);
# #   double              a_c = as<double>(a);
# #   vec    time_intervals_c = as<vec>(time_intervals);
# #   double t;
# # 
# #   vec    f_nm(Ns_c-1), g_nm(Ns_c-1), aux1(Ns_c-1), pow_gnm(Ns_c-1);
# #   double pi = 3.1415926535897932384626433832795028841971693993751058;
# # 
# #   vec        m = linspace<vec>(1, Ns_c-1, Ns_c-1);
# #   vec    n_all = linspace<vec>(1, Ns_c-1, Ns_c-1);
# #   double    k_0 = alpha_c*lambda_c*Ns_c;
# #   double    cte = 2*k_0*a_c/Ns_c;
# # 
# #   mat theory_timeseries1 = zeros<mat>(Ns_c, time_intervals_c.n_elem);
# #   g_nm    = 4*lambda_c*square(sin(0.5*pi*m/Ns_c));
# #   pow_gnm = pow(g_nm,-2);
# #   for(int n = 1; n<Ns_c; n++)
# #   {
# #     f_nm = sin(pi*m/Ns_c)%sin(pi*m*n/Ns_c);
# #     for(int i = 0; i<time_intervals_c.n_elem; i++)
# #     {
# #       t = time_intervals_c(i);
# #       theory_timeseries1(n-1, i) = as_scalar(sum(cte*f_nm%(pow_gnm - (t/g_nm + pow_gnm)%exp(-g_nm*t))));
# #     }  
# #   }
# #   
# #   for(int nn = 0; nn<aux1.n_elem; nn++) aux1(nn) = pow(-1., nn);
# #   f_nm = aux1%pow(cos(0.5*pi*m/Ns_c), 2);
# #   for(int i = 0; i<time_intervals_c.n_elem; i++)
# #   {
# #     t = time_intervals_c(i);
# #     theory_timeseries1(Ns_c-1, i) = as_scalar(sum(cte*f_nm%(pow(t,2)/2 - pow_gnm + (t/g_nm + pow_gnm)%exp(-g_nm*t))));
# #   }  
# #   return wrap(abs(theory_timeseries1));
# # "
# # fx.old <- cxxfunction(body = CA30Adenomas_AllClones.old.txt, sig = signature(a= "numeric", alpha= "numeric", lambda= "numeric", Ns= "integer",   time_intervals= "numeric"), plugin = "RcppArmadillo")
# 
# 
# 
# CA30Adenomas_AllClones.txt <- 
# "
#   using namespace arma;
#   int                Ns_c = as<int>(Ns);
#   double          alpha_c = as<double>(alpha);
#   double         lambda_c = as<double>(lambda);
#   double              a_c = as<double>(a);
#   vec    time_intervals_c = as<vec>(time_intervals);
#   double t;
# 
#   vec    f_nm(Ns_c-1), g_nm(Ns_c-1), aux1(Ns_c-1), pow_gnm(Ns_c-1);
#   double pi = 3.1415926535897932384626433832795028841971693993751058;
# 
#   vec        m = linspace<vec>(1, Ns_c-1, Ns_c-1);
#   vec    n_all = linspace<vec>(1, Ns_c-1, Ns_c-1);
#   double    k_0 = alpha_c*lambda_c*Ns_c;
#   double    cte = 2*k_0*a_c/Ns_c;
# 
#   mat theory_timeseries1 = zeros<mat>(Ns_c, time_intervals_c.n_elem);
#   g_nm    = 4*lambda_c*square(sin(0.5*pi*m/Ns_c));
#   pow_gnm = pow(g_nm,-2);
#   for(int n = 1; n<Ns_c; n++)
#   {
#     f_nm = sin(pi*m/Ns_c)%sin(pi*m*n/Ns_c);
#     for(int i = 0; i<time_intervals_c.n_elem; i++)
#     {
#       t = time_intervals_c(i);
#       theory_timeseries1(n-1, i) = as_scalar(sum(cte*f_nm%(t/g_nm - pow_gnm + pow_gnm%exp(-g_nm*t))));
#     }  
#   }
#   
#   for(int nn = 0; nn<aux1.n_elem; nn++) aux1(nn) = pow(-1., nn);
#   f_nm = aux1%pow(cos(0.5*pi*m/Ns_c), 2);
#   for(int i = 0; i<time_intervals_c.n_elem; i++)
#   {
#     t = time_intervals_c(i);
#     theory_timeseries1(Ns_c-1, i) = as_scalar(sum(cte*f_nm%(pow(t,2)/2 - t/g_nm + pow_gnm - pow_gnm%exp(-g_nm*t))));
#   }  
#   return wrap(abs(theory_timeseries1));
# "
# fx <- cxxfunction(body = CA30Adenomas_AllClones.txt, sig = signature(a= "numeric", alpha= "numeric", lambda= "numeric", Ns= "integer",   time_intervals= "numeric"), plugin = "RcppArmadillo")
# 
# CA30Adenomas_AllClones_minSize <- function(a, alpha, lambda, Ns, time.intervals, min.size = 4, full.size = 25)
# {
#   f_lambda_N.full <- fx(a, alpha, lambda, Ns, time.intervals)
#   min_N  <- round(min.size*Ns/full.size,0)
#   colSums(f_lambda_N.full[min_N:nrow(f_lambda_N.full),])
# }
# 
# CA30Adenomas_PartialClones <- cmpfun(CA30Adenomas_PartialClones)
# CA30Adenomas_FullClones    <- cmpfun(CA30Adenomas_FullClones)
# CA30Adenomas_AllClones     <- cmpfun(CA30Adenomas_AllClones)

# a              <- 6
# alpha          <- 1.1*10^(-4)
# lambda         <- 0.7
# Ns             <- 20
# time.intervals <- 1:400
# 
# system.time(uu1 <- CA30Adenomas_PartialClones(a, alpha, lambda, Ns, time.intervals))
# system.time(uu2 <- fx(a, alpha, lambda, Ns, time.intervals))
# system.time(uu3 <- CA30Adenomas_FullClones(a, alpha, lambda, Ns, time.intervals))
# system.time(uu4 <- CA30Adenomas_AllClones(a, alpha, lambda, Ns, time.intervals))
# print(max(uu1 - uu2[1:(Ns-1),]))
# print(max(uu3 - uu2[Ns,]))
# print(max(uu4 - colSums(uu2)))


