// Copyright (c) 2026 by Zijian Ye.
// Desc  : CPP code of MCMC for "Beyond Overall Diagnosis: A Longitudinal Multi-Domain Latent Structure Model for Domain-specific Parkinson’s Disease Progression".
// Note  : With the support of armspp, RcppArmadillo and RcppParallel.

#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include <limits>
#include <ctime>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "armspp.hpp"
#include <dqrng.h>
using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// [[Rcpp::depends(dqrng)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]

const double MIN_LOG_VALUE = -700.0; 
const double MIN_PROB = 1e-300;      
const double MAX_EXP_ARG = 700.0;    

inline double safe_log(double x)
{

    if (x <= MIN_PROB)
    {
        return MIN_LOG_VALUE;
    }
    return std::log(x);
}

inline double log1mexp(double x)
{

    if (x <= 0)
    {
        return MIN_LOG_VALUE;
    }
    if (x < 0.693)
    { 
        return std::log(-std::expm1(-x));
    }
    else
    {
        double exp_neg_x = std::exp(-x);
        if (exp_neg_x >= 1.0)
        {
            return MIN_LOG_VALUE;
        }
        return std::log1p(-exp_neg_x);
    }
}

// [[Rcpp::export]]
arma::vec f_rdirichlet(const arma::vec &alpha)
{

    size_t n = alpha.n_elem;
    arma::vec gamma_samples(n);

    for (size_t i = 0; i < n; ++i)
    {
        gamma_samples[i] = R::rgamma(alpha[i], 1.0);
    }

    double sum_gamma = arma::sum(gamma_samples);

    if (sum_gamma < 1e-100)
        return arma::ones<arma::vec>(n) / n;

    gamma_samples /= sum_gamma;

    return gamma_samples;
}

// [[Rcpp::export]]
double f_rwishart(int df, double sigma)
{

    double chi2_val = R::rchisq(df);

    return chi2_val * sigma;
}

// [[Rcpp::export]]
int f_rcategorical(const arma::vec &prob)
{

    double sum_prob = arma::sum(prob);

    double u = dqrng::dqrunif(1, 0.0, 1.0)[0] * sum_prob;

    double cum_prob = 0.0;

    for (size_t i = 0; i < prob.n_elem; ++i)
    {
        cum_prob += prob[i];
        if (u < cum_prob)
        {
            return i + 1;
        }
    }

    return prob.n_elem;
}

// [[Rcpp::export]]
arma::vec f_rmultinom(int n, const arma::vec &prob)
{

    arma::vec result(prob.n_elem, arma::fill::zeros);

    if (n == 0)
        return result;

    arma::vec cdf = arma::cumsum(prob);
    double max_prob = cdf(cdf.n_elem - 1);

    Rcpp::NumericVector u_vec = dqrng::dqrunif(n, 0.0, 1.0);

    for (int i = 0; i < n; ++i)
    {

        double u = u_vec[i] * max_prob;

        int index = 0;
        while (index < cdf.n_elem && u > cdf(index))
        {
            index++;
        }

        if (index >= prob.n_elem)
            index = prob.n_elem - 1;
        result[index]++;
    }

    return result;
}

// [[Rcpp::export]]
int f_rpoissonPositive(double lambda)
{

    if (lambda < 1e-6)
    {
        return 1;
    }

    int samp;
    do
    {
        samp = R::rpois(lambda);
    } while (samp == 0);

    return samp;
}

// [[Rcpp::export]]
arma::mat f_Mspline(int order, const arma::vec &x, const arma::vec &knots)
{

    int k1 = order;      
    int m = knots.n_elem;
    int n1 = m - 2 + k1; 
    arma::vec part1 = arma::ones<arma::vec>(k1) * knots(0);
    arma::vec part2 = knots.subvec(1, m - 2);
    arma::vec part3 = arma::ones<arma::vec>(k1) * knots(m - 1);
    arma::vec t1 = arma::join_cols(part1, part2, part3);

    arma::mat tem1 = arma::zeros<arma::mat>(n1 + k1 - 1, x.n_elem);
    for (int l = k1 - 1; l < n1; l++)
    {
        for (size_t i = 0; i < x.n_elem; i++)
        {
            if (x[i] >= t1[l] && x[i] < t1[l + 1])
            {
                tem1(l, i) = 1.0 / (t1[l + 1] - t1[l]);
            }
        }
    }

    arma::mat mbases;
    if (order == 1)
    {
        mbases = tem1;
    }
    else
    {
        mbases = tem1;
        for (int ii = 0; ii < order - 1; ii++)
        {
            arma::mat tem = arma::zeros<arma::mat>(n1 + k1 - 1 - (ii + 1), x.n_elem);
            for (int i = k1 - 1 - (ii + 1); i < n1; i++)
            {
                for (size_t j = 0; j < x.n_elem; j++)
                {
                    tem(i, j) = ((ii + 1) + 1) * ((x[j] - t1[i]) * mbases(i, j) + (t1[i + ii + 1 + 1] - x[j]) * mbases(i + 1, j)) /
                                ((t1[i + ii + 1 + 1] - t1[i]) * (ii + 1));
                }
            }
            mbases = tem;
        }
    }

    return mbases;
}

// [[Rcpp::export]]
arma::mat f_Ispline(int order, const arma::vec &x, const arma::vec &knots)
{

    int k = order + 1;    
    int m = knots.n_elem; 
    int n = m - 2 + k;    
    arma::vec t = arma::join_cols(arma::ones(k) * knots(0), knots.subvec(1, m - 2), arma::ones(k) * knots(m - 1));

    arma::mat yy1 = arma::zeros(n + k - 1, x.n_elem);

    for (int l = k - 1; l < n; l++)
    {
        for (size_t i = 0; i < x.n_elem; i++)
        {
            if (x(i) >= t(l) && x(i) < t(l + 1))
            {
                yy1(l, i) = 1.0 / (t(l + 1) - t(l));
            }
        }
    }

    arma::mat yytem1 = yy1;
    for (int ii = 1; ii <= order; ii++)
    {
        arma::mat yytem2 = arma::zeros(n + k - 1 - ii, x.n_elem);
        for (int i = k - ii - 1; i < n; i++)
        {
            for (size_t j = 0; j < x.n_elem; j++)
            {
                yytem2(i, j) = (ii + 1) * ((x(j) - t(i)) * yytem1(i, j) + (t(i + ii + 1) - x(j)) * yytem1(i + 1, j)) /
                               (t(i + ii + 1) - t(i)) / ii;
            }
        }
        yytem1 = yytem2;
    }

    arma::vec index = arma::zeros(x.n_elem);
    for (size_t i = 0; i < x.n_elem; i++)
    {
        index(i) = arma::sum(t <= x(i));
    }

    arma::mat yy = arma::zeros(n - 1, x.n_elem);

    if (order == 1)
    {
        for (int i = 1; i < n; i++)
        {
            for (size_t j = 0; j < x.n_elem; j++)
            {
                yy(i - 1, j) = ((i + 1) < index(j) - order + 1) +
                               ((i + 1) == index(j)) * (t(i + order + 1) - t(i)) * yytem1(i, j) / (order + 1);
            }
        }
    }
    else
    {
        for (size_t j = 0; j < x.n_elem; j++)
        {
            for (int i = 1; i < n; i++)
            {
                if ((i + 1) < (index(j) - order + 1))
                {
                    yy(i - 1, j) = 1;
                }
                else if ((i + 1) <= index(j) && (i + 1) >= (index(j) - order + 1))
                {
                    yy(i - 1, j) = arma::dot(t.subvec(i + order + 1, index(j) + order) - t.subvec(i, index(j) - 1),
                                             yytem1.col(j).subvec(i, index(j) - 1)) /
                                   (order + 1);
                }
            }
        }
    }

    return yy;
}

// [[Rcpp::export]]
arma::vec f_Ispline_x(int order, double x, const arma::vec &knots)
{

    int k = order + 1;    
    int m = knots.n_elem; 
    int n = m - 2 + k;    
    arma::vec t = arma::join_cols(arma::ones(k) * knots(0), knots.subvec(1, m - 2), arma::ones(k) * knots(m - 1));

    arma::vec yy1 = arma::zeros(n + k - 1);

    for (int l = k - 1; l < n; l++)
    {
        if (x >= t(l) && x < t(l + 1))
        {
            yy1(l) = 1.0 / (t(l + 1) - t(l));
        }
    }

    arma::vec yytem1 = yy1;
    for (int ii = 1; ii <= order; ii++)
    {
        arma::vec yytem2 = arma::zeros(n + k - 1 - ii);
        for (int i = k - ii - 1; i < n; i++)
        {
            yytem2(i) = (ii + 1) * ((x - t(i)) * yytem1(i) + (t(i + ii + 1) - x) * yytem1(i + 1)) /
                        (t(i + ii + 1) - t(i)) / ii;
        }
        yytem1 = yytem2;
    }

    double index = arma::sum(t <= x);

    arma::vec yy = arma::zeros(n - 1);

    if (order == 1)
    {
        for (int i = 1; i < n; i++)
        {
            yy(i - 1) = ((i + 1) < index - order + 1) +
                        ((i + 1) == index) * (t(i + order + 1) - t(i)) * yytem1(i) / (order + 1);
        }
    }
    else
    {
        for (int i = 1; i < n; i++)
        {
            if ((i + 1) < (index - order + 1))
            {
                yy(i - 1) = 1;
            }
            else if ((i + 1) <= index && (i + 1) >= (index - order + 1))
            {
                yy(i - 1) = arma::dot(t.subvec(i + order + 1, index + order) - t.subvec(i, index - 1),
                                      yytem1.subvec(i, index - 1)) /
                            (order + 1);
            }
        }
    }

    return yy;
}

// [[Rcpp::export]]
arma::mat f_Ispline_Inter(int N, int G, int order, int Lknots, const arma::mat &x, const arma::vec &knots)
{

    arma::mat bis(Lknots, N * G, arma::fill::zeros);

    for (int g = 0; g < G; ++g)
    {
        arma::mat bis2 = f_Ispline(order, x.col(g), knots);
        for (int i = 0; i < N; ++i)
        {
            bis.col(i + N * g) = bis2.col(i);
        }
    }

    return bis;
}

// [[Rcpp::export]]
arma::mat f_R_na(const arma::mat &InterR)
{

    arma::mat R2 = InterR;

    for (arma::uword i = 0; i < R2.n_rows; ++i)
    {
        for (arma::uword j = 0; j < R2.n_cols; ++j)
        {
            if (R2(i, j) == 10000)
            {
                R2(i, j) = 0;
            }
        }
    }
    return R2;
}

// [[Rcpp::export]]
arma::mat f_FRig(int N, int G, int order, const arma::vec &ri, const arma::mat &gammaL, const arma::mat &beta, const arma::mat &X, const arma::mat &tig, const arma::vec &knots, const arma::vec &Ui)
{

    arma::mat FRig(arma::max(ri) + 1, G * N, arma::fill::zeros);

    arma::mat bis(gammaL.n_elem / G, N * (max(ri) + 1));
    double stash_exp;
    double LambdatL;
    double LambdatR;
    double FInterL;
    double FInterR;

    for (int i = 0; i < N; ++i)
    {
        for (int t = 0; t < ri(i); ++t)
        {
            bis.col(i + N * t) = f_Ispline_x(order, tig(t, i), knots);
        }
    }

    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            stash_exp = exp(arma::dot(beta.col(g).t(), X.col(i + N * g)) + Ui(i));

            FRig(0, g + G * i) = 1 - exp(-arma::dot(gammaL.col(g), bis.col(i)) * stash_exp);

            FRig(ri(i), g + G * i) = 1 - (1 - exp(-arma::dot(gammaL.col(g), bis.col(i + N * (ri(i) - 1))) * stash_exp));

            for (int t = 1; t < ri(i); ++t)
            {

                LambdatL = arma::dot(gammaL.col(g), bis.col(i + N * (t - 1)));
                LambdatR = arma::dot(gammaL.col(g), bis.col(i + N * t));

                FInterL = 1 - exp(-LambdatL * stash_exp);
                FInterR = 1 - exp(-LambdatR * stash_exp);

                FRig(t, g + G * i) = FInterR - FInterL;
            }
        }
    }

    return FRig;
}

// [[Rcpp::export]]
arma::mat f_logFRig(int N, int G, int order, const arma::vec &ri, const arma::mat &gammaL, const arma::mat &beta, const arma::mat &X, const arma::mat &tig, const arma::vec &knots, const arma::vec &Ui)
{

    arma::mat logFRig(arma::max(ri) + 1, G * N, arma::fill::zeros);

    arma::mat bis(gammaL.n_elem / G, N * (max(ri) + 1));
    double stash_exp;
    double LambdatL;
    double LambdatR;
    double FInterL;
    double FInterR;

    for (int i = 0; i < N; ++i)
    {
        for (int t = 0; t < ri(i); ++t)
        {
            bis.col(i + N * t) = f_Ispline_x(order, tig(t, i), knots);
        }
    }

    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            stash_exp = exp(arma::dot(beta.col(g).t(), X.col(i + N * g)) + Ui(i));
            double lambda0 = arma::dot(gammaL.col(g), bis.col(i)) * stash_exp;

            logFRig(0, g + G * i) = log1mexp(lambda0);

            logFRig(ri(i), g + G * i) = -arma::dot(gammaL.col(g), bis.col(i + N * (ri(i) - 1))) * stash_exp;

            for (int t = 1; t < ri(i); ++t)
            {

                LambdatL = arma::dot(gammaL.col(g), bis.col(i + N * (t - 1)));
                LambdatR = arma::dot(gammaL.col(g), bis.col(i + N * t));

                FInterL = 1 - exp(-LambdatL * stash_exp);
                FInterR = 1 - exp(-LambdatR * stash_exp);

                double diff = FInterR - FInterL;
                logFRig(t, g + G * i) = safe_log(std::max(diff, MIN_PROB));
            }
        }
    }

    return logFRig;
}

// [[Rcpp::export]]
arma::mat f_PRigt(int N, int J, int G, const arma::vec &s, const arma::vec &ri, const arma::mat &Y, const arma::mat &Theta, const arma::mat &logFRig)
{

    int max_ri = ri.max();
    int max_ri_1 = (max_ri + 1);
    arma::mat PRigt(G, (max_ri + 1) * N, arma::fill::zeros);
    double deno1 = 0.0;
    double deno2 = 0.0;
    arma::mat deno_log(max_ri + 1, N, arma::fill::zeros);
    arma::uvec indices;
    arma::vec logFRig_col;
    arma::vec deno_log_col;
    double m;
    double deno3;

    for (int g = 0; g < G; ++g)
    {
        indices = arma::find(s == g + 1);
        for (int i = 0; i < N; ++i)
        {
            deno1 = 0;
            deno2 = 0;
            deno_log.zeros();

            for (int t = 0; t < ri[i] + 1; ++t)
            {
                for (auto j : indices)
                {
                    for (int r = 0; r < t; ++r)
                    {
                        if (Y(j, r + max_ri * i) != 0)
                        {
                            deno_log(t, i) += safe_log(Theta(Y(j, r + max_ri * i) - 1, 0 + 2 * j));
                        }
                    }
                    for (int r = t; r < ri[i]; ++r)
                    {
                        if (Y(j, r + max_ri * i) != 0)
                        {
                            deno_log(t, i) += safe_log(Theta(Y(j, r + max_ri * i) - 1, 1 + 2 * j));
                        }
                    }
                }
            }

            logFRig_col = logFRig.col(g + G * i).head(ri[i] + 1);
            deno_log_col = deno_log.col(i).head(ri[i] + 1);
            m = arma::max(deno_log_col + logFRig_col);

            for (int t = 0; t < ri[i] + 1; ++t)
            {
                deno1 += std::exp(logFRig(t, g + G * i) + deno_log(t, i) - m);
            }

            deno2 = m + safe_log(deno1);

            for (int t = 0; t < ri[i] + 1; ++t)
            {
                deno3 = deno_log(t, i) + logFRig(t, g + G * i);
                PRigt(g, t + max_ri_1 * i) = std::exp(deno3 - deno2);
            }
        }
    }

    return PRigt;
}

// [[Rcpp::export]]
arma::mat f_Rig(int N, int G, const arma::vec &ri, const arma::mat &PRigt)
{

    int max_ri = ri.max();
    arma::mat R(N, G, arma::fill::zeros);

    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            arma::vec prob = (PRigt.row(g).subvec(i * (max_ri + 1), i * (max_ri + 1) + ri[i])).t();

            R(i, g) = f_rcategorical(prob);
        }
    }

    return R;
}

// [[Rcpp::export]]
arma::mat f_InterL(int N, int G, const arma::vec &ri, const arma::mat &R, const arma::mat &tig)
{

    arma::mat InterL(N, G, arma::fill::zeros);

    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            int Rig = R(i, g);
            if (Rig == 1)
            {
                InterL(i, g) = 0;
            }
            else
            {
                InterL(i, g) = tig(Rig - 2, i);
            }
        }
    }

    return InterL;
}

// [[Rcpp::export]]
arma::mat f_InterR(int N, int G, const arma::vec &ri, const arma::mat &R, const arma::mat &tig)
{

    arma::mat InterR(N, G, arma::fill::zeros);

    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            int Rig = R(i, g);
            if (Rig < ri(i) + 1)
            {
                InterR(i, g) = tig(Rig - 1, i);
            }
            else
            {
                InterR(i, g) = 10000;
            }
        }
    }

    return InterR;
}

// [[Rcpp::export]]
arma::mat f_status(int N, int G, const arma::vec &ri, const arma::mat &R, const arma::mat &tig)
{

    arma::mat status(N, G, arma::fill::zeros);

    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            int Rig = R(i, g);
            if (Rig == 1)
            {
                status(i, g) = 0;
            }
            else if (Rig == ri(i) + 1)
            {
                status(i, g) = 2;
            }
            else
            {
                status(i, g) = 1;
            }
        }
    }

    return status;
}

// [[Rcpp::export]]
arma::vec f_XI(const arma::mat &sg)
{

    arma::vec alpha3 = arma::sum(sg, 1) + 1.0;
    arma::vec XI = f_rdirichlet(alpha3);

    return XI;
}

// [[Rcpp::export]]
arma::mat f_Psg(int N, int J, int G, const arma::vec &ri, const arma::mat &R, const arma::mat &Theta, const arma::mat &Y, const arma::vec &XI)
{

    int rimax = ri.max();
    arma::mat Psg(G, J, arma::fill::zeros);
    arma::vec deno_log(G, arma::fill::zeros);
    int ri_max_i;
    int R_ig;
    int Y_val;
    arma::vec log_XI;
    double m;
    double deno1;
    double deno2;

    for (int j = 0; j < J; ++j)
    {
        deno_log.zeros();

        for (int i = 0; i < N; ++i)
        {
            ri_max_i = rimax * i;
            for (int g = 0; g < G; ++g)
            {
                R_ig = R(i, g);
                for (int t = 0; t < R_ig - 1; ++t)
                {
                    Y_val = Y(j, t + ri_max_i);
                    if (Y_val != 0)
                    {
                        deno_log(g) += safe_log(Theta(Y_val - 1, 0 + 2 * j));
                    }
                }
                for (int t = R_ig - 1; t < ri[i]; ++t)
                {
                    Y_val = Y(j, t + ri_max_i);
                    if (Y_val != 0)
                    {
                        deno_log(g) += safe_log(Theta(Y_val - 1, 1 + 2 * j));
                    }
                }
            }
        }

        log_XI = arma::log(XI);
        m = arma::max(deno_log + log_XI);
        deno1 = 0.0;

        for (int g = 0; g < G; ++g)
        {
            deno1 += std::exp(log_XI(g) + deno_log(g) - m);
        }

        deno2 = m + safe_log(deno1);

        for (int g = 0; g < G; ++g)
        {
            Psg(g, j) = std::exp(log_XI(g) + deno_log(g) - deno2);
        }
    }

    return Psg;
}

struct PsgWorker : public Worker
{

    const int N;
    const int J;
    const int G;
    const arma::vec &ri;
    const arma::mat &R;
    const arma::mat &Theta;
    const arma::mat &Y;
    const arma::vec &XI;

    arma::mat &Psg;

    PsgWorker(int N, int J, int G, const arma::vec &ri, const arma::mat &R, const arma::mat &Theta, const arma::mat &Y, const arma::vec &XI, arma::mat &Psg)
        : N(N), J(J), G(G), ri(ri), R(R), Theta(Theta), Y(Y), XI(XI), Psg(Psg) {}

    void operator()(std::size_t j_start, std::size_t j_end)
    {
        int rimax = ri.max();
        arma::vec deno_log(G, arma::fill::zeros);
        int ri_max_i;
        int R_ig;
        int Y_val;
        arma::vec log_XI;
        double m;
        double deno1;
        double deno2;

        for (std::size_t j = j_start; j < j_end; ++j)
        {
            deno_log.zeros();

            for (int i = 0; i < N; ++i)
            {
                ri_max_i = rimax * i;
                for (int g = 0; g < G; ++g)
                {
                    R_ig = R(i, g);
                    for (int t = 0; t < R_ig - 1; ++t)
                    {
                        Y_val = Y(j, t + ri_max_i);
                        if (Y_val != 0)
                        {
                            deno_log(g) += safe_log(Theta(Y_val - 1, 0 + 2 * j));
                        }
                    }
                    for (int t = R_ig - 1; t < ri[i]; ++t)
                    {
                        Y_val = Y(j, t + ri_max_i);
                        if (Y_val != 0)
                        {
                            deno_log(g) += safe_log(Theta(Y_val - 1, 1 + 2 * j));
                        }
                    }
                }
            }

            log_XI = arma::log(XI);
            m = arma::max(deno_log + log_XI);
            deno1 = 0.0;

            for (int g = 0; g < G; ++g)
            {
                deno1 += std::exp(log_XI(g) + deno_log(g) - m);
            }

            deno2 = m + safe_log(deno1);

            for (int g = 0; g < G; ++g)
            {
                Psg(g, j) = std::exp(log_XI(g) + deno_log(g) - deno2);
            }
        }
    }
};

// [[Rcpp::export]]
arma::vec f_sj(int J, int G, int ir, int BI, const arma::mat &Psg)
{

    arma::vec s(J);

    for (int j = 0; j < J; ++j)
    {
        s[j] = f_rcategorical(Psg.col(j));
    }

    if (ir + 1 < 0.1 * BI)
    {
        if ((ir + 1) % static_cast<int>(0.01 * BI) == 0)
        {
            arma::vec prob(G, arma::fill::ones);
            prob /= G;

            for (int j = 0; j < J; ++j)
            {
                s[j] = f_rcategorical(prob);
            }
        }
    }

    arma::vec unique_s = arma::unique(s);
    if (unique_s.n_elem < G)
    {
        arma::vec prob(G, arma::fill::ones);
        prob /= G;

        for (int j = 0; j < J; ++j)
        {
            s[j] = f_rcategorical(prob);
        }
    }

    return s;
}

// [[Rcpp::export]]
arma::mat f_sg(int J, int G, const arma::vec &s)
{

    arma::mat sg(G, J, arma::fill::zeros);

    for (int j = 0; j < J; ++j)
    {
        sg(s[j] - 1, j) = 1;
    }

    return sg;
}

// [[Rcpp::export]]
arma::mat f_Theta(int N, int J, const arma::vec &ri, const arma::vec &dj, const arma::vec &s, const arma::mat &R, const arma::mat &Yijtc)
{

    int max_dj = dj.max();
    int max_ri = ri.max();
    arma::mat Theta(max_dj, 2 * J, arma::fill::zeros);
    int g;
    int R_ig;
    int idx_dj;
    int max_ri_j;

    for (int j = 0; j < J; ++j)
    {
        g = s[j] - 1;
        max_ri_j = max_ri * j;
        arma::vec v1(dj[j], arma::fill::ones);
        arma::vec v2(dj[j], arma::fill::ones);
        for (int dj_count = 0; dj_count < dj[j]; ++dj_count)
        {
            for (int i = 0; i < N; ++i)
            {
                R_ig = R(i, g);
                idx_dj = dj_count + max_dj * i;
                for (int t = 0; t < R_ig - 1; ++t)
                {
                    v1[dj_count] += Yijtc(t + max_ri_j, idx_dj);
                }
                for (int t = R_ig - 1; t < ri[i]; ++t)
                {
                    v2[dj_count] += Yijtc(t + max_ri_j, idx_dj);
                }
            }
        }

        Theta.col(0 + 2 * j).head(dj[j]) = f_rdirichlet(v1);
        Theta.col(1 + 2 * j).head(dj[j]) = f_rdirichlet(v2);
    }

    return Theta;
}

// [[Rcpp::export]]
arma::mat f_label_switching_Theta(int N, int J, int G, const arma::vec &s, arma::mat Theta)
{

    for (int g = 0; g < G; ++g)
    {
        arma::uvec sgroup = arma::find(s == (g + 1));
        int indexg = sgroup.n_elem;

        if (indexg > 0)
        {
            arma::vec prob = arma::ones<arma::vec>(indexg) / indexg;
            int indexj = f_rcategorical(prob) - 1;

            if (Theta(0, 2 * sgroup[indexj]) < Theta(0, 2 * sgroup[indexj] + 1))
            {
                for (int j = 0; j < indexg; ++j)
                {
                    int col0 = 2 * sgroup[j];
                    int col1 = col0 + 1;
                    Theta.cols(col0, col1).swap_cols(0, 1);
                }
            }
        }
    }

    return Theta;
}

// [[Rcpp::export]]
arma::mat f_LambdatL(int N, int G, const arma::mat &gammaL, const arma::mat &bisL)
{

    arma::mat LambdatL(N, G, arma::fill::zeros);

    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            LambdatL(i, g) = arma::dot(gammaL.col(g), bisL.col(i + N * g));
        }
    }

    return LambdatL;
}

// [[Rcpp::export]]
arma::mat f_LambdatR(int N, int G, const arma::mat &gammaL, const arma::mat &bisR)
{

    arma::mat LambdatR(N, G, arma::fill::zeros);

    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            LambdatR(i, g) = arma::dot(gammaL.col(g), bisR.col(i + N * g));
        }
    }

    return LambdatR;
}

// [[Rcpp::export]]
arma::mat f_z(int N, int G, const arma::mat &status, const arma::mat &LambdatR, const arma::mat &X, const arma::mat &beta, const arma::vec &Ui)
{

    arma::mat z(N, G, arma::fill::zeros);
    double templam1;

    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            if (status(i, g) == 0)
            {
                templam1 = LambdatR(i, g) * std::exp(arma::dot(beta.col(g).t(), X.col(i + N * g)) + Ui(i));
                z(i, g) = f_rpoissonPositive(templam1);
            }
        }
    }

    return z;
}

// [[Rcpp::export]]
arma::mat f_zz(int N, int G, int Lknots, const arma::mat &status, const arma::mat &LambdatR, const arma::mat &X, const arma::mat &beta, const arma::mat &gammaL, const arma::mat &bisR, const arma::mat &z)
{

    arma::mat zz(Lknots, N * G, arma::fill::zeros);

    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            if (status(i, g) == 0)
            {
                zz.col(i + N * g) = f_rmultinom(z(i, g), gammaL.col(g) % bisR.col(i + N * g));
            }
        }
    }

    return zz;
}

// [[Rcpp::export]]
arma::mat f_w(int N, int G, const arma::mat &status, const arma::mat &LambdatR, const arma::mat &LambdatL, const arma::mat &X, const arma::mat &beta, const arma::vec &Ui)
{

    arma::mat w(N, G, arma::fill::zeros);
    double templam1;

    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            if (status(i, g) == 1)
            {
                templam1 = (LambdatR(i, g) - LambdatL(i, g)) * std::exp(arma::dot(beta.col(g).t(), X.col(i + N * g)) + Ui(i));
                w(i, g) = f_rpoissonPositive(templam1);
            }
        }
    }

    return w;
}

// [[Rcpp::export]]
arma::mat f_ww(int N, int G, int Lknots, const arma::mat &status, const arma::mat &LambdatR, const arma::mat &LambdatL, const arma::mat &X, const arma::mat &beta, const arma::mat &gammaL, const arma::mat &bisR, const arma::mat &bisL, const arma::mat &w)
{

    arma::mat ww(Lknots, N * G, arma::fill::zeros);

    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            if (status(i, g) == 1)
            {
                ww.col(i + N * g) = f_rmultinom(w(i, g), gammaL.col(g) % (bisR.col(i + N * g) - bisL.col(i + N * g)));
            }
        }
    }

    return ww;
}

// [[Rcpp::export]]
arma::vec f_te1(const arma::vec &z, const arma::vec &w, const arma::vec status)
{

    arma::vec te1(z.n_elem, arma::fill::zeros);

    for (arma::uword i = 0; i < z.n_elem; ++i)
    {
        if (status[i] == 0)
        {
            te1[i] = z[i];
        }
        else if (status[i] == 1)
        {
            te1[i] = w[i];
        }
    }

    return te1;
}

// [[Rcpp::export]]
arma::vec f_te2(const arma::vec &LambdatR, const arma::vec &LambdatL, const arma::vec &status)
{

    arma::vec te2(LambdatR.n_elem, arma::fill::zeros);

    for (arma::uword i = 0; i < LambdatR.n_elem; ++i)
    {
        if (status[i] == 0 || status[i] == 1)
        {
            te2[i] = LambdatR[i];
        }
        else if (status[i] == 2)
        {
            te2[i] = LambdatR[i] + LambdatL[i];
        }
    }

    return te2;
}

// [[Rcpp::export]]
arma::vec f_beta_C(arma::vec beta, const arma::mat &xx, const arma::vec &te1, const arma::vec &te2, double sig0, double coef_range, const arma::vec &Ui)
{

    struct BetaLogPdf
    {
        arma::mat xx;
        arma::vec beta;
        arma::vec te1;
        arma::vec te2;
        arma::vec Ui;
        double sig0;
        int j;

        BetaLogPdf(const arma::mat &xx_, int j_, arma::vec beta_, const arma::vec &te1_, const arma::vec &te2_, double sig0_, const arma::vec &Ui_)
            : xx(xx_), beta(beta_), te1(te1_), te2(te2_), Ui(Ui_), sig0(sig0_), j(j_) {}

        double operator()(double x) const
        {
            arma::vec beta_copy = beta;
            beta_copy(j) = x;
            double result = (arma::sum((xx.col(j) * x + Ui) % te1) -
                             arma::sum(arma::exp(xx * beta_copy + Ui) % te2) -
                             x * x / (2 * sig0 * sig0));
            return result;
        }
    };

    int nInitial = 10;
    arma::vec Initial = -coef_range + linspace(1, 10, 10) * (coef_range + coef_range) / (nInitial + 1);
    double convexity = 0;

    int maxPoints = max(2 * nInitial + 1, 100);
    bool metropolis = false;

    int p = beta.n_elem;

    for (int j = 0; j < p; ++j)
    {

        BetaLogPdf logPdf(xx, j, beta, te1, te2, sig0, Ui);

        armspp::ARMS<double, BetaLogPdf, arma::vec::iterator> sampler(
            logPdf, -coef_range, coef_range, convexity, Initial.begin(), nInitial, maxPoints, metropolis, beta[j]);

        beta[j] = sampler();
    }

    return beta;
}

// [[Rcpp::export]]
arma::mat f_Beta(int N, int G, int B, arma::mat beta, const arma::mat &X, double sig0, double coef_range, const arma::mat &z, const arma::mat &w, const arma::mat &status, const arma::mat &LambdatR, const arma::mat &LambdatL, const arma::vec &Ui)
{

    for (int g = 0; g < G; ++g)
    {
        arma::mat xx(N, B);
        for (int i = 0; i < N; ++i)
        {
            xx.row(i) = X.col(i + N * g).t();
        }
        arma::vec te1 = f_te1(z.col(g), w.col(g), status.col(g));
        arma::vec te2 = f_te2(LambdatR.col(g), LambdatL.col(g), status.col(g));
        beta.col(g) = f_beta_C(beta.col(g), xx, te1, te2, sig0, coef_range, Ui);
    }

    return beta;
}

// [[Rcpp::export]]
arma::mat f_gammaL(int N, int G, int Lknots, const arma::mat &zz, const arma::mat &ww, const arma::mat &status, const arma::mat &bisR, const arma::mat &bisL, const arma::mat &X, const arma::mat &beta, const arma::vec &eta, const arma::vec &Ui)
{

    arma::mat gammaL(Lknots, G);

    for (int g = 0; g < G; ++g)
    {
        for (int l = 0; l < Lknots; ++l)
        {
            double tempa = 1.0;
            double tempb = eta(g);

            for (arma::uword i = 0; i < N; ++i)
            {
                if (status(i, g) == 0 && bisR(l, i + N * g) > 0)
                {
                    tempa += zz(l, i + N * g);
                }
                else if (status(i, g) == 1 && (bisR(l, i + N * g) - bisL(l, i + N * g)) > 0)
                {
                    tempa += ww(l, i + N * g);
                }
            }

            for (arma::uword i = 0; i < N; ++i)
            {
                if (status(i, g) == 0)
                {
                    tempb += bisR(l, i + N * g) * std::exp(arma::dot(X.col(i + N * g), beta.col(g)) + Ui(i));
                }
                else if (status(i, g) == 1)
                {
                    tempb += bisR(l, i + N * g) * std::exp(arma::dot(X.col(i + N * g), beta.col(g)) + Ui(i));
                }
                else if (status(i, g) == 2)
                {
                    tempb += bisL(l, i + N * g) * std::exp(arma::dot(X.col(i + N * g), beta.col(g)) + Ui(i));
                }
            }

            gammaL(l, g) = R::rgamma(tempa, 1.0 / tempb);
        }
    }

    return gammaL;
}

// [[Rcpp::export]]
arma::vec f_eta(int G, int Lknots, arma::mat gammaL, double a_eta, double b_eta)
{

    arma::vec eta(G);

    for (int g = 0; g < G; ++g)
    {
        eta(g) = R::rgamma(a_eta + Lknots, 1.0 / (b_eta + arma::sum(gammaL.col(g))));
    }

    return eta;
}

// [[Rcpp::export]]
double f_SigmaG(int N, double sigma_alpha, double sigma_beta, const arma::vec &Ui)
{

    double sigma_alpha1 = sigma_alpha + 0.5 * N;
    double sigma_beta1 = sigma_beta + 0.5 * arma::dot(Ui, Ui);

    double SigmaG = 1.0 / R::rgamma(sigma_alpha1, 1 / sigma_beta1);

    return SigmaG;
}

// [[Rcpp::export]]
double f_Ui_C(double Ui, arma::mat beta, const arma::vec &xx, const arma::vec &te1, const arma::vec &te2, double sig0, double coef_range)
{

    struct UiLogPdf
    {
        arma::vec xx;
        arma::mat beta;
        arma::vec te1;
        arma::vec te2;
        double sig0;

        UiLogPdf(const arma::vec &xx_, const arma::mat &beta_, const arma::vec &te1_, const arma::vec &te2_, double sig0_)
            : xx(xx_), beta(beta_), te1(te1_), te2(te2_), sig0(sig0_) {}

        double operator()(double x) const
        {
            double result = (arma::sum(x * te1) -
                             arma::sum(std::exp(x) * (arma::exp((xx.t() * beta).t()) % te2)) -
                             x * x / (2 * sig0));
            return result;
        }
    };

    int nInitial = 10;
    arma::vec Initial = -coef_range + linspace(1, 10, 10) * (coef_range + coef_range) / (nInitial + 1);
    double convexity = 0;

    int maxPoints = max(2 * nInitial + 1, 100);
    bool metropolis = false;

    UiLogPdf logPdf(xx, beta, te1, te2, sig0);

    armspp::ARMS<double, UiLogPdf, arma::vec::iterator> sampler(
        logPdf, -coef_range, coef_range, convexity, Initial.begin(), nInitial, maxPoints, metropolis, Ui);

    Ui = sampler();

    return Ui;
}

// [[Rcpp::export]]
arma::vec f_Ui(arma::vec Ui, int N, int G, int B, arma::mat beta, const arma::mat &X, double sig0, double coef_range, const arma::mat &z, const arma::mat &w, const arma::mat status, const arma::mat &LambdatR, const arma::mat &LambdatL)
{

    for (int i = 0; i < N; ++i)
    {
        arma::vec xx(B);

        xx = X.col(i + N * 0);
        arma::vec te1 = f_te1(z.row(i).t(), w.row(i).t(), status.row(i).t());
        arma::vec te2 = f_te2(LambdatR.row(i).t(), LambdatL.row(i).t(), status.row(i).t());

        Ui[i] = f_Ui_C(Ui[i], beta, xx, te1, te2, sig0, coef_range);
    }

    return Ui;
}

// [[Rcpp::export]]
arma::vec f_knots(const arma::vec &InterL, const arma::vec &R2, int length = 10)
{

    double min_val = std::min(arma::min(InterL), arma::min(R2));
    double max_val = std::max(arma::max(InterL), arma::max(R2));

    arma::vec knots = arma::linspace(min_val, max_val, length);

    return knots;
}

// [[Rcpp::export]]
double f_log_Likelihood(int N, int J, int G, int order, int Lknots, const arma::vec &ri, const arma::mat &Y, const arma::vec &s, const arma::mat &Theta, const arma::mat &gammaL, const arma::vec &knots, const arma::mat &beta, const arma::mat &X, const arma::mat &tig, const arma::mat &R, const arma::vec &Ui)
{

    double L = 0.0;
    int max_ri = max(ri);
    arma::mat logFRig = f_logFRig(N, G, order, ri, gammaL, beta, X, tig, knots, Ui);

    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            int R_ig = R(i, g);
            for (int j = 0; j < J; ++j)
            {
                if (s[j] == g + 1)
                {
                    for (int r = 0; r < R_ig - 1; ++r)
                    {
                        if (Y(j, r + max_ri * i) != 0)
                        {
                            L += safe_log(Theta(Y(j, r + max_ri * i) - 1, 0 + 2 * j));
                        }
                    }
                    for (int r = R_ig - 1; r < ri[i]; ++r)
                    {
                        if (Y(j, r + max_ri * i) != 0)
                        {
                            L += safe_log(Theta(Y(j, r + max_ri * i) - 1, 1 + 2 * j));
                        }
                    }
                }
            }
            L += logFRig(R_ig - 1, g + G * i);
        }
    }

    return L;
}

// [[Rcpp::export]]
arma::vec f_log_py(int N, int J, int G, int order, int Lknots, const arma::vec &ri, const arma::mat &Y, const arma::vec &s, const arma::mat &Theta, const arma::mat &gammaL, const arma::vec &knots, const arma::mat &beta, const arma::mat &X, const arma::mat &tig, const arma::mat &R, const arma::vec &Ui)
{

    arma::vec log_py(N, arma::fill::zeros);
    int max_ri = max(ri);
    arma::mat logFRig = f_logFRig(N, G, order, ri, gammaL, beta, X, tig, knots, Ui);

    for (int i = 0; i < N; ++i)
    {
        int max_ri_i = max_ri * i;
        for (int g = 0; g < G; ++g)
        {
            int R_ig_1 = R(i, g) - 1;
            for (int j = 0; j < J; ++j)
            {
                int j2 = 2 * j;
                int j2_1 = 1 + 2 * j;
                if (s[j] == g + 1)
                {
                    for (int r = 0; r < R_ig_1; ++r)
                    {
                        if (Y(j, r + max_ri_i) != 0)
                        {
                            log_py(i) += safe_log(Theta(Y(j, r + max_ri_i) - 1, j2));
                        }
                    }
                    for (int r = R_ig_1; r < ri[i]; ++r)
                    {
                        if (Y(j, r + max_ri_i) != 0)
                        {
                            log_py(i) += safe_log(Theta(Y(j, r + max_ri_i) - 1, j2_1));
                        }
                    }
                }
            }
            log_py(i) += logFRig(R_ig_1, g + G * i);
        }
    }

    return log_py;
}

// [[Rcpp::export]]
double f_AIC(int N, int U, const arma::mat &Ulog_py, int pd)
{

    double AIC = 0.0;

    AIC = -2.0 * (1.0 / U) * (accu(Ulog_py)) + 2.0 * pd;

    return AIC;
}

// [[Rcpp::export]]
double f_BIC(int N, int U, int pd, double LKLH_mean)
{

    double BIC = 0.0;

    BIC = pd * std::log(N) - 2.0 * LKLH_mean;

    return BIC;
}

// [[Rcpp::export]]
double f_WAIC(int N, int U, const arma::mat &Ulog_p_yx)
{
    
    double WAIC = 0.0;
    double lppd = 0.0;  
    double pwaic = 0.0; 

    for (int i = 0; i < N; i++)
    {
        double alpha = max(Ulog_p_yx.col(i));

        double col_var = var(Ulog_p_yx.col(i), 0);
        if (std::isnan(col_var) || std::isinf(col_var))
        {
            col_var = 0.0;
        }

        double log_mean_exp = -std::log(U) + alpha + safe_log(sum(exp(Ulog_p_yx.col(i) - alpha)));

        lppd += log_mean_exp;
        pwaic += col_var;
    }

    WAIC = -2.0 * (lppd - pwaic);

    return WAIC;
}

//[[Rcpp::export]]
int f_pd(int N, int J, int G, int Lknotes, int B, const arma::vec &dj)
{

    int pd = 2 * sum(dj - 1) + G * Lknotes + B * G + 1 + (N * G + J);

    return pd;
}

// [[Rcpp::export]]
int f_mode_sorted(int size, const arma::vec &TOC_mode)
{

    int mode = TOC_mode(0);
    int current = 1;
    int mostOccured = 1;

    for (int x = 1; x < size; ++x)
    {
        if (TOC_mode(x) == TOC_mode(x - 1))
        {
            current++;
        }
        else
        {
            if (current > mostOccured)
            {
                mostOccured = current;
                mode = TOC_mode(x - 1);
            }
            current = 1;
        }
    }

    if (current > mostOccured)
    {
        mostOccured = current;
        mode = TOC_mode(size - 1);
    }

    return mode;
}

// [[Rcpp::export]]
List f_mcmc(int IR, int BI, int N, int J, int G, int B, int order, int Lknots, arma::vec ri, arma::vec dj, arma::mat Y, arma::vec s, arma::vec knots, arma::mat Theta, arma::mat Yijtc, arma::mat beta, arma::mat gammaL, arma::vec eta, arma::mat X, arma::mat tig, arma::vec Ui, double sig0, double coef_range_beta, double coef_range_ui, double a_eta, double b_eta, double SigmaG, int sigma_alpha, double sigma_beta, int ds)
{

    arma::mat logFRig(max(ri + 1), G * N, arma::fill::ones);
    arma::mat PRigt(G, max(ri + 1) * N, arma::fill::zeros);
    arma::mat R(N, G, arma::fill::zeros);
    arma::mat InterR(N, G, arma::fill::zeros);
    arma::mat InterR2(N, G, arma::fill::zeros);
    arma::mat InterL(N, G, arma::fill::zeros);
    arma::mat status(N, G, arma::fill::zeros);
    arma::mat sg(G, J, arma::fill::zeros);
    arma::vec XI(G, arma::fill::zeros);
    arma::mat Psg(G, J, arma::fill::zeros);
    arma::mat z(N, G, arma::fill::zeros);
    arma::mat w(N, G, arma::fill::zeros);
    arma::mat zz(Lknots, N * G, arma::fill::zeros);
    arma::mat ww(Lknots, N * G, arma::fill::zeros);
    arma::mat LambdatR(N, G, arma::fill::zeros);
    arma::mat LambdatL(N, G, arma::fill::zeros);
    arma::mat bisR(Lknots, N * G, arma::fill::zeros);
    arma::mat bisL(Lknots, N * G, arma::fill::zeros);

    arma::mat s_e(J, IR, arma::fill::zeros);
    arma::mat Theta_e(max(dj), 2 * J * IR, arma::fill::zeros);
    arma::mat R_e(N, G * IR, arma::fill::zeros);
    arma::mat beta_e(B, G * IR, arma::fill::zeros);
    arma::mat gammaL_e(Lknots, G * IR, arma::fill::zeros);
    arma::mat PRigt_e(G, max(ri + 1) * N, arma::fill::zeros);
    arma::mat Ui_e(N, IR, arma::fill::zeros);
    arma::vec SigmaG_e(IR, arma::fill::zeros);
    arma::mat eta_e(G, IR, arma::fill::zeros);

    arma::vec count_censoring(G, arma::fill::zeros);
    arma::vec count_first(G, arma::fill::zeros);
    arma::mat Ulog_py(IR - BI + 1, N, arma::fill::zeros);

    int pd = 0;
    double BIC = 0.0;
    double AIC = 0.0;
    double WAIC = 0.0;
    double LogLikelihood = 0.0;

    logFRig = f_logFRig(N, G, order, ri, gammaL, beta, X, tig, knots, Ui);
    PRigt = f_PRigt(N, J, G, s, ri, Y, Theta, logFRig);
    R = f_Rig(N, G, ri, PRigt);

    InterL = f_InterL(N, G, ri, R, tig);
    InterR = f_InterR(N, G, ri, R, tig);
    status = f_status(N, G, ri, R, tig);
    InterR2 = f_R_na(InterR);

    bisL = f_Ispline_Inter(N, G, order, Lknots, InterL, knots);
    bisR = f_Ispline_Inter(N, G, order, Lknots, InterR2, knots);
    LambdatL = f_LambdatL(N, G, gammaL, bisL);
    LambdatR = f_LambdatR(N, G, gammaL, bisR);

    for (int ir = 0; ir < IR; ++ir)
    {

        Theta = f_Theta(N, J, ri, dj, s, R, Yijtc);
        Theta = f_label_switching_Theta(N, J, G, s, Theta);

        sg = f_sg(J, G, s);
        XI = f_XI(sg);
        Psg = f_Psg(N, J, G, ri, R, Theta, Y, XI);
        s = f_sj(J, G, ir, BI, Psg);

        logFRig = f_logFRig(N, G, order, ri, gammaL, beta, X, tig, knots, Ui);
        PRigt = f_PRigt(N, J, G, s, ri, Y, Theta, logFRig);
        R = f_Rig(N, G, ri, PRigt);

        InterL = f_InterL(N, G, ri, R, tig);
        InterR = f_InterR(N, G, ri, R, tig);
        status = f_status(N, G, ri, R, tig);
        InterR2 = f_R_na(InterR);

        bisL = f_Ispline_Inter(N, G, order, Lknots, InterL, knots);
        bisR = f_Ispline_Inter(N, G, order, Lknots, InterR2, knots);

        z = f_z(N, G, status, LambdatR, X, beta, Ui);
        zz = f_zz(N, G, Lknots, status, LambdatR, X, beta, gammaL, bisR, z);
        w = f_w(N, G, status, LambdatR, LambdatL, X, beta, Ui);
        ww = f_ww(N, G, Lknots, status, LambdatR, LambdatL, X, beta, gammaL, bisR, bisL, w);

        beta = f_Beta(N, G, B, beta, X, sig0, coef_range_beta, z, w, status, LambdatR, LambdatL, Ui);

        gammaL = f_gammaL(N, G, Lknots, zz, ww, status, bisR, bisL, X, beta, eta, Ui);

        LambdatL = f_LambdatL(N, G, gammaL, bisL);
        LambdatR = f_LambdatR(N, G, gammaL, bisR);

        eta = f_eta(G, Lknots, gammaL, a_eta, b_eta);

        SigmaG = f_SigmaG(N, sigma_alpha, sigma_beta, Ui);
        Ui = f_Ui(Ui, N, G, B, beta, X, SigmaG, coef_range_ui, z, w, status, LambdatR, LambdatL);

        s_e.col(ir) = s;
        eta_e.col(ir) = eta;
        SigmaG_e(ir) = SigmaG;
        Ui_e.col(ir) = Ui;
        for (int k = 0; k < 2; ++k)
        {
            for (int j = 0; j < J; ++j)
            {
                Theta_e.col(k + 2 * j + 2 * J * ir).head(dj[j]) = Theta.col(k + 2 * j).head(dj[j]);
            }
        }
        for (int g = 0; g < G; ++g)
        {
            R_e.col(g + G * ir) = R.col(g);
        }
        for (int g = 0; g < G; ++g)
        {
            gammaL_e.col(g + G * ir) = gammaL.col(g);
        }
        for (int g = 0; g < G; ++g)
        {
            beta_e.col(g + G * ir) = beta.col(g);
        }
    }

    PRigt_e /= (IR - BI + 1);

    arma::mat mean_Theta(max(dj), 2 * J, arma::fill::zeros);
    for (int djcount = 0; djcount < max(dj); ++djcount)
    {
        for (int k = 0; k < 2; ++k)
        {
            for (int j = 0; j < J; ++j)
            {
                for (int ir = BI - 1; ir < IR; ++ir)
                {
                    mean_Theta(djcount, k + 2 * j) += Theta_e(djcount, k + 2 * j + 2 * J * ir);
                }
            }
        }
    }
    mean_Theta /= (IR - BI + 1);

    arma::vec mean_Ui(N, arma::fill::zeros);
    for (int i = 0; i < N; ++i)
    {
        for (int ir = BI - 1; ir < IR; ++ir)
        {
            mean_Ui(i) += Ui_e(i, ir);
        }
    }
    mean_Ui /= (IR - BI + 1);

    arma::mat mean_gammaL(Lknots, G, arma::fill::zeros);
    for (int l = 0; l < Lknots; ++l)
    {
        for (int g = 0; g < G; ++g)
        {
            for (int ir = BI - 1; ir < IR; ++ir)
            {
                mean_gammaL(l, g) += gammaL_e(l, g + G * ir);
            }
        }
    }
    mean_gammaL /= (IR - BI + 1);

    double mean_SigmaG = 0.0;
    for (int ir = BI - 1; ir < IR; ++ir)
    {
        mean_SigmaG += SigmaG_e(ir);
    }
    mean_SigmaG /= (IR - BI + 1);

    double median_SigmaG = arma::median(SigmaG_e.subvec(BI, IR - 1));

    arma::vec mean_eta(G, arma::fill::zeros);
    for (int ir = BI - 1; ir < IR; ++ir)
    {
        mean_eta += eta_e.col(ir);
    }
    mean_eta /= (IR - BI + 1);

    arma::mat mode_R(N, G, arma::fill::zeros);
    arma::vec TOC_mode_R(IR - BI + 1, arma::fill::zeros);
    for (int i = 0; i < N; ++i)
    {
        for (int g = 0; g < G; ++g)
        {
            for (int ir = BI - 1; ir < IR; ++ir)
            {
                TOC_mode_R(ir - BI + 1) = R_e(i, g + G * ir);
            }
            TOC_mode_R = sort(TOC_mode_R, "ascend");
            mode_R(i, g) = f_mode_sorted(IR - BI + 1, TOC_mode_R);
        }
    }

    arma::vec TOC_mode_S(IR - BI + 1, arma::fill::zeros);
    arma::vec mode_S(J, arma::fill::zeros);
    for (int j = 0; j < J; ++j)
    {
        TOC_mode_S = sort(s_e.row(j).subvec(BI - 1, IR - 1), "ascend").t();
        mode_S(j) = f_mode_sorted(IR - BI + 1, TOC_mode_S);
    }

    for (int g = 0; g < G; g++)
    {
        for (int i = 0; i < N; i++)
        {
            if (mode_R.col(g)[i] == 1)
                count_first[g] += 1;
            if (mode_R.col(g)[i] == (ri[i] + 1))
                count_censoring[g] += 1;
        }
    }

    arma::mat lower_ci(B, G, arma::fill::zeros);
    arma::mat upper_ci(B, G, arma::fill::zeros);
    arma::mat mean_beta(B, G, arma::fill::zeros);
    arma::mat sd_beta(B, G, arma::fill::zeros);
    arma::cube beta_samples(IR - BI + 1, B, G);

    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            for (int ir = BI - 1; ir < IR; ++ir)
            {
                beta_samples(ir - (BI - 1), b, g) = beta_e(b, g + G * ir);
            }
        }
    }

    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {

            arma::vec samples = beta_samples.slice(g).col(b);
            samples = arma::sort(samples);

            mean_beta(b, g) = arma::mean(samples);
            sd_beta(b, g) = arma::stddev(samples);
            int lower_index = static_cast<int>(0.025 * samples.n_elem);
            int upper_index = static_cast<int>(0.975 * samples.n_elem);
            lower_ci(b, g) = samples(lower_index);
            upper_ci(b, g) = samples(upper_index);
        }
    }
   
    LogLikelihood = f_log_Likelihood(N, J, G, order, Lknots, ri, Y, mode_S, mean_Theta, mean_gammaL, knots, mean_beta, X, tig, mode_R, mean_Ui);

    return List::create(
        Rcpp::Named("s_e") = s_e,
        Rcpp::Named("Theta_e") = Theta_e,
        Rcpp::Named("gammaL_e") = gammaL_e,
        Rcpp::Named("beta_e") = beta_e,
        Rcpp::Named("PRigt_e") = PRigt_e,
        Rcpp::Named("SigmaG_e") = SigmaG_e,
        Rcpp::Named("mode_S") = mode_S,
        Rcpp::Named("mode_R") = mode_R,
        Rcpp::Named("mean_Theta") = mean_Theta,
        Rcpp::Named("mean_eta") = mean_eta,
        Rcpp::Named("mean_beta") = mean_beta,
        Rcpp::Named("mean_gammaL") = mean_gammaL,
        Rcpp::Named("mean_Ui") = mean_Ui,
        Rcpp::Named("mean_SigmaG") = mean_SigmaG,
        Rcpp::Named("median_SigmaG") = median_SigmaG,
        Rcpp::Named("BIC") = BIC,
        Rcpp::Named("AIC") = AIC,
        Rcpp::Named("WAIC") = WAIC,
        Rcpp::Named("LogLikelihood") = LogLikelihood);
}

// [[Rcpp::export]]
void f_beta_show(int B, int G, int IR, int BI, arma::mat beta_e)
{

    arma::mat lower_ci(B, G, arma::fill::zeros);
    arma::mat upper_ci(B, G, arma::fill::zeros);
    arma::mat mean_beta(B, G, arma::fill::zeros);
    arma::mat sd_beta(B, G, arma::fill::zeros);
    arma::cube beta_samples(IR - BI + 1, B, G);

    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            for (int ir = BI - 1; ir < IR; ++ir)
            {
                beta_samples(ir - (BI - 1), b, g) = beta_e(b, g + G * ir);
            }
        }
    }

    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {

            arma::vec samples = beta_samples.slice(g).col(b);
            samples = arma::sort(samples);

            mean_beta(b, g) = arma::mean(samples);
            sd_beta(b, g) = arma::stddev(samples);
            int lower_index = static_cast<int>(0.025 * samples.n_elem);
            int upper_index = static_cast<int>(0.975 * samples.n_elem);
            lower_ci(b, g) = samples(lower_index);
            upper_ci(b, g) = samples(upper_index);
        }
    }

    cout << "Beta(coef): " << endl;
    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            std::cout << "Beta[" << b + 1 << "][" << g + 1 << "] - "
                      << "Mean: " << std::setw(5) << std::right << std::setfill(' ') << mean_beta(b, g)
                      << ", SD: " << std::setw(5) << std::right << std::setfill(' ') << sd_beta(b, g)
                      << ", CI: [" << std::setw(5) << std::right << std::setfill(' ') << lower_ci(b, g) << ", "
                      << std::setw(5) << std::right << std::setfill(' ') << upper_ci(b, g) << "]";
            if (lower_ci(b, g) > 0 || upper_ci(b, g) < 0)
            {
                std::cout << " *";
            }
            std::cout << std::endl;
        }
    }
}

// [[Rcpp::export]]
arma::cube f_gammaL_CI(int B, int G, int Lknots, int IR, int BI, arma::mat gammaL_e)
{

    arma::mat lower_ci(Lknots, G, arma::fill::zeros);
    arma::mat upper_ci(Lknots, G, arma::fill::zeros);
    arma::mat mean_gammaL(Lknots, G, arma::fill::zeros);
    arma::mat sd_gammaL(Lknots, G, arma::fill::zeros);
    arma::cube gammaL_samples(IR - BI + 1, Lknots, G);

    for (int l = 0; l < Lknots; ++l)
    {
        for (int g = 0; g < G; ++g)
        {
            for (int ir = BI - 1; ir < IR; ++ir)
            {
                gammaL_samples(ir - (BI - 1), l, g) = gammaL_e(l, g + G * ir);
            }
        }
    }

    for (int l = 0; l < Lknots; ++l)
    {
        for (int g = 0; g < G; ++g)
        {

            arma::vec samples = gammaL_samples.slice(g).col(l);
            samples = arma::sort(samples);

            mean_gammaL(l, g) = arma::mean(samples);
            sd_gammaL(l, g) = arma::stddev(samples);
            int lower_index = static_cast<int>(0.025 * samples.n_elem);
            int upper_index = static_cast<int>(0.975 * samples.n_elem);
            lower_ci(l, g) = samples(lower_index);
            upper_ci(l, g) = samples(upper_index);
        }
    }

    arma::cube ci_cube(lower_ci.n_rows, lower_ci.n_cols, 2);
    ci_cube.slice(0) = lower_ci;
    ci_cube.slice(1) = upper_ci;

    return ci_cube;
}

